(() => {
    const SCORE_DESCRIPTIONS = {
        AlphaMissense: "AlphaMissense score (0-1) trained on protein context and evolutionary information.",
        CADD_phred: "CADD phred score. Values above 20 indicate variants in the top 1% of predicted deleteriousness.",
        REVEL_score: "REVEL ensemble score. Values above 0.5 suggest damaging effect.",
        MetaRNN_score: "MetaRNN recurrent neural network score where values above 0.5 indicate likely damaging.",
        MetaSVM_score: "MetaSVM meta-predictor integrating multiple functional annotations.",
        MetaLR_score: "MetaLR logistic regression score combining annotation features.",
        BayesDel_addAF_score: "BayesDel score that incorporates allele frequency information.",
        BayesDel_noAF_score: "BayesDel score without allele frequency information.",
        ClinPred_score: "ClinPred pathogenicity score for missense variants.",
        MVP_score: "MVP score for missense pathogenicity prediction.",
        MutPred_score: "MutPred score estimating molecular impact of amino acid substitutions.",
        MPC_score: "MPC score quantifying regional missense constraint.",
        PrimateAI_score: "PrimateAI score using deep learning on primate variation.",
        DEOGEN2_score: "DEOGEN2 score integrating sequence and interaction features.",
        "M-CAP_score": "M-CAP score tuned for rare missense pathogenicity classification.",
        VEST4_score: "VEST4 functional impact score.",
        DANN_score: "DANN deep neural network score for variant deleteriousness.",
        LINSIGHT: "LINSIGHT score estimating noncoding fitness consequences.",
        FATHMM_score: "FATHMM score based on hidden Markov models.",
        PROVEAN_score: "PROVEAN functional effect score where lower values suggest deleteriousness.",
        SIFT_score: "SIFT score where values below 0.05 are predicted damaging.",
        Polyphen2_HDIV_score: "PolyPhen-2 HumDiv score; higher values indicate probable damage.",
        Polyphen2_HVAR_score: "PolyPhen-2 HumVar score for Mendelian-like interpretation.",
        LRT_score: "Likelihood ratio test score for selective constraint.",
        MutationTaster_score: "MutationTaster pathogenicity confidence score.",
        MutationAssessor_score: "MutationAssessor functional impact score.",
        "GERP++_RS": "GERP++ rejected substitution score reflecting evolutionary constraint.",
        phyloP100way_vertebrate: "phyloP conservation score across 100 vertebrates.",
        phastCons100way_vertebrate: "phastCons probability of conserved element.",
        SiPhy_29way_logOdds: "SiPhy log odds score from 29-way alignments.",
        integrated_fitCons_score: "fitCons integrated score for functional consequence likelihood.",
    };

    const chartRegistry = {};
    let currentSingleResult = null;
    let currentBatchPayload = null;
    let currentBatchPage = 1;
    let batchSortState = { key: "variant", asc: true };
    let batchFilterState = "all";
    let batchGeneFilter = "";
    let currentGeneResults = [];

    function safeValue(value) {
        return value === null || value === undefined || value === "" ? "--" : value;
    }

    function asNumber(value) {
        if (value === null || value === undefined || value === "" || value === "--") return null;
        const num = Number(value);
        return Number.isFinite(num) ? num : null;
    }

    function normalizeScore(name, value) {
        const num = asNumber(value);
        if (num === null) return null;
        if (name === "CADD_phred") return Math.min(1, Math.max(0, num / 40));
        if (name === "SIFT_score") return Math.min(1, Math.max(0, 1 - num));
        if (name === "PROVEAN_score") return Math.min(1, Math.max(0, (Math.abs(num) / 12)));
        return Math.min(1, Math.max(0, num));
    }

    // Cut-offs follow the dbNSFP / tool documentation. dir: "ge" = damaging when value >= cut, "le" = when <= cut.
    const TOOL_CUTOFFS = {
        CADD_phred: ["ge", 20], SIFT_score: ["le", 0.05], PROVEAN_score: ["le", -2.5], FATHMM_score: ["le", -1.5],
        LRT_score: ["le", 0.001], BayesDel_addAF_score: ["ge", 0.0692655], BayesDel_noAF_score: ["ge", -0.0570105],
        "M-CAP_score": ["ge", 0.025], MetaSVM_score: ["ge", 0], MetaLR_score: ["ge", 0.5], MetaRNN_score: ["ge", 0.5],
        REVEL_score: ["ge", 0.5], ClinPred_score: ["ge", 0.5], PrimateAI_score: ["ge", 0.803], DEOGEN2_score: ["ge", 0.5],
        Polyphen2_HDIV_score: ["ge", 0.453], Polyphen2_HVAR_score: ["ge", 0.447], MutationAssessor_score: ["ge", 1.935],
        MPC_score: ["ge", 2], MVP_score: ["ge", 0.7], AlphaMissense: ["ge", 0.564], MutPred_score: ["ge", 0.5],
        VEST4_score: ["ge", 0.5], DANN_score: ["ge", 0.5],
    };

    function predictionForScore(name, value) {
        const num = asNumber(value);
        if (num === null) return "--";
        const rule = TOOL_CUTOFFS[name];
        if (!rule) return "(no standard cut-off)";
        const damaging = rule[0] === "ge" ? num >= rule[1] : num <= rule[1];
        return `${damaging ? "Damaging" : "Tolerated"} (${rule[0] === "ge" ? "≥" : "≤"} ${rule[1]})`;
    }

    function destroyChart(name) {
        if (chartRegistry[name]) {
            chartRegistry[name].destroy();
            delete chartRegistry[name];
        }
    }

    function drawChart(name, canvasId, config) {
        const canvas = document.getElementById(canvasId);
        if (!canvas || !window.Chart) return;
        destroyChart(name);
        chartRegistry[name] = new Chart(canvas, config);
    }

    async function postForm(endpoint, formData) {
        const response = await fetch(endpoint, { method: "POST", body: formData });
        const payload = await response.json();
        if (!response.ok || payload.status === "error") {
            throw new Error(payload.message || "Request failed");
        }
        return payload;
    }

    async function postJson(endpoint, payload) {
        const response = await fetch(endpoint, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        });
        const data = await response.json();
        if (!response.ok || data.status === "error") {
            throw new Error(data.message || "Request failed");
        }
        return data;
    }

    function showInlineMessage(target, message, kind = "info") {
        if (!target) return;
        const classes = kind === "error"
            ? "border-red-200 bg-red-50 text-red-700"
            : "border-teal-200 bg-teal-50 text-teal-800";
        target.innerHTML = `<div class="rounded-lg border ${classes} px-3 py-2 text-sm">${message}</div>`;
    }

    async function singleVariantSearch(chr, pos, ref, alt) {
        const preview = document.getElementById("single-result-preview");
        if (preview) {
            preview.innerHTML = "<div class='text-sm text-neutral'>Searching variant...</div>";
        }

        const formData = new FormData();
        formData.set("chrom", chr);
        formData.set("pos", pos);
        formData.set("ref", ref);
        formData.set("alt", alt);

        try {
            const data = await postForm("/api/v1/variant", formData);
            currentSingleResult = data;
            renderResults(data);

            const url = new URL(window.location.href);
            url.searchParams.set("chr", data.variant.chr);
            url.searchParams.set("pos", data.variant.pos);
            url.searchParams.set("ref", data.variant.ref);
            url.searchParams.set("alt", data.variant.alt);
            window.history.replaceState({}, "", url.toString());
            return data;
        } catch (error) {
            if (preview) {
                showInlineMessage(preview, error.message === "Variant not found" ? "Not among the precomputed GRCh38 missense alleles on this server (no score; absence is not evidence of benignity)." : escapeHtml(error.message), "error");
            }
            throw error;
        }
    }

    function updateBadge(element, text) {
        if (!element) return;
        element.textContent = text || "--";
        element.classList.remove("badge-pathogenic", "badge-benign", "badge-vus");
        const value = String(text || "").toLowerCase();
        if (value.includes("pathogenic") || value.includes("damaging")) {
            element.classList.add("badge-pathogenic");
        } else if (value.includes("benign") || value.includes("tolerated")) {
            element.classList.add("badge-benign");
        } else {
            element.classList.add("badge-vus");
        }
    }

    function escapeHtml(value) {
        return String(value ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
    }

    function fmtScore(value, digits = 4) {
        const num = asNumber(value);
        return num === null || value === null || value === undefined ? "--" : num.toFixed(digits);
    }

    function matchBadge(value) {
        if (value === true) return '<span class="rounded-full bg-green-100 px-2 py-0.5 text-xs font-semibold text-green-800">yes</span>';
        if (value === false) return '<span class="rounded-full bg-amber-100 px-2 py-0.5 text-xs font-semibold text-amber-800">no — different protein/isoform</span>';
        return '<span class="text-xs text-neutral">n/a</span>';
    }

    function renderTranscripts(tx) {
        const body = document.getElementById("txTableBody");
        const badge = document.getElementById("txPerAllele");
        if (!body) return;
        if (!tx) {
            body.innerHTML = '<tr><td colspan="4" class="px-3 py-2 text-neutral">Transcript context needs the annotation database.</td></tr>';
            return;
        }
        if (badge) {
            badge.textContent = tx.annotatemissense_score_is_per_allele
                ? `1 score for this allele · ${tx.annotation_rows} database row${tx.annotation_rows === 1 ? "" : "s"}`
                : `${tx.distinct_annotatemissense_scores} distinct scores across rows`;
        }
        const d = tx.dbnsfp || {};
        const am = tx.alphamissense || {};
        const positions = (d.aa_positions || []).map((p) => String(p).replace(/\.0$/, "")).join("/");
        const rows = [];
        rows.push(`<tr><td class="px-3 py-2 font-semibold">dbNSFP v5.1 (challenge input)</td>
            <td class="px-3 py-2 font-mono text-xs">${escapeHtml((d.transcripts || []).join(", ") || "--")}<br><span class="text-neutral">${escapeHtml((d.genes || []).join(", "))}</span></td>
            <td class="px-3 py-2 font-mono text-xs">${escapeHtml(d.ref_aa || "?")}${positions ? escapeHtml(positions) : "&gt;"}${escapeHtml(d.alt_aa || "?")}</td>
            <td class="px-3 py-2 text-xs text-neutral">reference</td></tr>`);
        rows.push(`<tr><td class="px-3 py-2 font-semibold">AlphaMissense</td>
            <td class="px-3 py-2 font-mono text-xs">${escapeHtml(am.transcript || "--")}<br><span class="text-neutral">${escapeHtml(am.uniprot || "")}</span></td>
            <td class="px-3 py-2 font-mono text-xs">${escapeHtml(am.protein_variant || "--")}</td>
            <td class="px-3 py-2">${matchBadge(am.matches_dbnsfp_amino_acids)}</td></tr>`);
        (tx.esm1v || []).forEach((e, i) => {
            rows.push(`<tr><td class="px-3 py-2 font-semibold">ESM1v${tx.esm1v.length > 1 ? ` (row ${i + 1})` : ""}</td>
                <td class="px-3 py-2 font-mono text-xs">${escapeHtml(e.transcript || "--")}<br><span class="text-neutral">${escapeHtml(e.protein || "")}</span></td>
                <td class="px-3 py-2 font-mono text-xs">${escapeHtml(e.hgvs_p || "--")}</td>
                <td class="px-3 py-2">${matchBadge(e.matches_dbnsfp_amino_acids)}</td></tr>`);
        });
        if (!(tx.esm1v || []).length) {
            rows.push('<tr><td class="px-3 py-2 font-semibold">ESM1v</td><td colspan="3" class="px-3 py-2 text-xs text-neutral">No ESM1v record joined for this allele.</td></tr>');
        }
        rows.push(`<tr><td class="px-3 py-2 font-semibold">ANNOVAR refGene</td>
            <td colspan="3" class="px-3 py-2 font-mono text-xs">${(tx.annovar_refgene_aachange || []).map(escapeHtml).join("<br>") || "--"}</td></tr>`);
        body.innerHTML = rows.join("");
    }

    function renderResults(data) {
        currentSingleResult = data;
        const am = data.annotatemissense || {};

        const preview = document.getElementById("single-result-preview");
        if (preview) {
            const v = data.variant;
            const cls = am.prediction === "D" ? "border-red-200 bg-red-50" : am.prediction === "T" ? "border-green-200 bg-green-50" : "border-slate-200 bg-slate-50";
            preview.innerHTML = `
                <div class="rounded-lg border ${cls} p-3 text-sm">
                    <p class="font-mono font-semibold">${escapeHtml(v.chr)}:${escapeHtml(v.pos)} ${escapeHtml(v.ref)}&gt;${escapeHtml(v.alt)} · ${escapeHtml(data.gene?.name || "")}</p>
                    <p class="mt-1">AnnotateMissense score <strong>${fmtScore(am.score)}</strong> · ${escapeHtml(am.label)} · ClinVar: ${escapeHtml(data.clinical?.CLNSIG || "--")}</p>
                    <a class="mt-2 inline-block font-semibold text-primary underline" href="/results?chr=${encodeURIComponent(v.chr)}&pos=${encodeURIComponent(v.pos)}&ref=${encodeURIComponent(v.ref)}&alt=${encodeURIComponent(v.alt)}">Open full results →</a>
                </div>
            `;
        }

        const variantTitle = document.getElementById("variantHeaderTitle");
        if (!variantTitle) return;

        document.getElementById("resultsFallback")?.classList.add("hidden");
        ["section-model", "section-info", "section-transcripts", "section-scores", "section-population", "section-conservation", "section-compare"].forEach((id) => {
            document.getElementById(id)?.classList.remove("hidden");
        });

        const variant = data.variant || {};
        const gene = data.gene || {};
        const clinical = data.clinical || {};
        const esm = data.esm1v || {};
        const population = data.population || {};
        const conservation = data.conservation || {};
        const scores = data.pathogenicity_scores || {};

        variantTitle.textContent = `${safeValue(variant.chr)}:${safeValue(variant.pos)} ${safeValue(variant.ref)} → ${safeValue(variant.alt)}`;
        document.getElementById("variantHeaderMeta").textContent = `GRCh38 · ${safeValue(gene.name)} · ${safeValue(gene.function)}`;
        const scoreOnly = data.provenance && data.provenance.annotation_rows === null;
        updateBadge(document.getElementById("clinvarBadge"), clinical.CLNSIG ? `ClinVar: ${clinical.CLNSIG}` : scoreOnly ? "ClinVar: not loaded here" : "Not in ClinVar");
        if (scoreOnly) {
            ["section-info", "section-transcripts", "section-scores", "section-population", "section-conservation"].forEach((id) => document.getElementById(id)?.classList.add("hidden"));
            document.getElementById("esmChart")?.parentElement?.classList.add("hidden");
            document.getElementById("download-json")?.classList.add("hidden");
        }

        const scoreNum = asNumber(am.score);
        document.getElementById("amScoreText").textContent = fmtScore(am.score, 6);
        document.getElementById("amPredText").textContent = am.prediction ? `${am.prediction} (${am.label})` : "--";
        document.getElementById("amConfText").textContent = fmtScore(am.confidence, 6);
        document.getElementById("amRowsText").textContent = safeValue(data.provenance?.annotation_rows);
        document.getElementById("amGaugeValue").textContent = fmtScore(am.score, 4);
        document.getElementById("amGaugeLabel").textContent = am.label || "--";
        document.getElementById("combinedScoreText").textContent = fmtScore(esm.combined_score, 3);

        const gaugeColor = scoreNum === null ? "#94a3b8" : scoreNum >= 0.5 ? "#dc2626" : "#15803d";
        drawChart("amGauge", "amGauge", {
            type: "doughnut",
            data: {
                labels: ["AnnotateMissense score", "Remaining"],
                datasets: [{
                    data: [Math.max(0, Math.min(1, scoreNum || 0)), 1 - Math.max(0, Math.min(1, scoreNum || 0))],
                    backgroundColor: [gaugeColor, "#e2e8f0"],
                    borderWidth: 0,
                }],
            },
            options: { cutout: "72%", plugins: { legend: { display: false }, tooltip: { enabled: false } } },
        });

        drawChart("esm", "esmChart", {
            type: "bar",
            data: {
                labels: ["Model 1", "Model 2", "Model 3", "Model 4", "Model 5"],
                datasets: [{ data: (esm.models || []).map((v) => asNumber(v) || 0), backgroundColor: "#0f766e" }],
            },
            options: { plugins: { legend: { display: false } }, scales: { y: { title: { display: true, text: "log-likelihood ratio" } } } },
        });

        renderTranscripts(data.transcripts);

        document.getElementById("infoRsid").textContent = safeValue(variant.rsid);
        document.getElementById("infoGene").textContent = safeValue(gene.name);
        document.getElementById("infoCytoband").textContent = safeValue(gene.cytoband);
        const hg19 = variant.hg19 || {};
        document.getElementById("infoHg19").textContent = hg19.pos ? `chr${String(hg19.chr || "").replace(/^chr/, "")}:${hg19.pos}` : "--";
        document.getElementById("infoAa").textContent = safeValue(gene.aa_change);
        document.getElementById("infoFunction").textContent = safeValue(gene.function);
        document.getElementById("infoDomain").textContent = safeValue(gene.interpro_domain);
        const ensemblLink = document.getElementById("infoEnsembl");
        const ensemblGene = safeValue(gene.ensembl_gene);
        ensemblLink.textContent = ensemblGene;
        if (ensemblGene !== "--") {
            ensemblLink.href = `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(ensemblGene)}`;
        }

        const clinSig = document.getElementById("clinSigText");
        updateBadge(clinSig, clinical.CLNSIG || "Unspecified");
        document.getElementById("clinDisease").textContent = clinical.CLNDN ? String(clinical.CLNDN).replace(/_/g, " ").split("|").join("; ") : "--";
        document.getElementById("clinReview").textContent = clinical.CLNREVSTAT ? String(clinical.CLNREVSTAT).replace(/_/g, " ").replace(/ ,/g, ",") : "--";

        const reviewText = String(clinical.CLNREVSTAT || "").toLowerCase();
        let stars = 1;
        if (reviewText.includes("practice_guideline")) stars = 4;
        else if (reviewText.includes("multiple_submitters")) stars = 3;
        else if (reviewText.includes("single_submitter")) stars = 2;
        document.getElementById("clinStars").textContent = "★".repeat(stars) + "☆".repeat(4 - stars);

        const clinvarLink = document.getElementById("clinvarLink");
        if (variant.rsid) {
            clinvarLink.href = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(variant.rsid)}`;
        } else {
            clinvarLink.href = "https://www.ncbi.nlm.nih.gov/clinvar/";
        }

        const scoreEntries = Object.entries(scores).filter(([k]) => !k.endsWith("_class") && !k.endsWith("_pred"));
        const scoreChartKeys = ["AlphaMissense", "CADD_phred", "REVEL_score", "MetaRNN_score", "BayesDel_addAF_score"];
        const chartVals = scoreChartKeys.map((name) => normalizeScore(name, scores[name]));
        drawChart("scoreComparison", "scoreComparisonChart", {
            type: "bar",
            data: {
                labels: scoreChartKeys,
                datasets: [{
                    data: chartVals.map((v) => v ?? 0),
                    backgroundColor: chartVals.map((v) => (v !== null && v >= 0.5 ? "#dc2626" : "#16a34a")),
                }],
            },
            options: {
                scales: { y: { min: 0, max: 1 } },
                plugins: { legend: { display: false } },
            },
        });

        const tableBody = document.getElementById("scoreTableBody");
        tableBody.innerHTML = "";
        scoreEntries.forEach(([name, rawValue], idx) => {
            const normalized = normalizeScore(name, rawValue);
            const width = normalized !== null ? Math.round(normalized * 100) : 0;
            const prediction = predictionForScore(name, rawValue);
            const tr = document.createElement("tr");
            tr.className = "score-main-row cursor-pointer";
            tr.dataset.scoreName = name;
            tr.innerHTML = `
                <td class="px-3 py-2 font-medium">${name}</td>
                <td class="px-3 py-2">${safeValue(rawValue)}</td>
                <td class="px-3 py-2">${prediction}</td>
                <td class="px-3 py-2">
                    <div class="score-bar-track">
                        <div class="score-bar-fill" style="width:${width}%;background:${prediction.startsWith("Damaging") ? "#dc2626" : prediction.startsWith("Tolerated") ? "#16a34a" : "#94a3b8"}"></div>
                    </div>
                </td>
            `;
            const detail = document.createElement("tr");
            detail.className = "score-row-detail hidden";
            detail.dataset.detailFor = name;
            detail.innerHTML = `<td class="px-3 py-2 text-xs text-neutral" colspan="4">${SCORE_DESCRIPTIONS[name] || "No description available."}</td>`;
            tr.addEventListener("click", () => {
                detail.classList.toggle("hidden");
            });
            tableBody.appendChild(tr);
            tableBody.appendChild(detail);
            if (idx > 200) return;
        });

        const popVals = [
            asNumber(population.gnomad30_AF_afr) || 0,
            asNumber(population.gnomad30_AF_amr) || 0,
            asNumber(population.gnomad30_AF_eas) || 0,
            asNumber(population.gnomad30_AF_nfe) || 0,
            asNumber(population.gnomad30_AF_sas) || 0,
            asNumber(population.gnomad30_AF_fin) || 0,
            asNumber(population.gnomad30_AF) || 0,
        ];
        drawChart("population", "populationChart", {
            type: "bar",
            data: {
                labels: ["AFR", "AMR", "EAS", "NFE", "SAS", "FIN", "Overall"],
                datasets: [{ data: popVals, backgroundColor: "#0f766e" }],
            },
            options: { plugins: { legend: { display: false } } },
        });

        document.getElementById("popNote").textContent = popVals.every((v) => v === 0)
            ? "Not observed in gnomAD"
            : "Observed in gnomAD populations";
        document.getElementById("exacText").textContent = safeValue(population.ExAC_ALL);
        document.getElementById("popmaxText").textContent = safeValue(population.AF_popmax);

        drawChart("conservation", "conservationChart", {
            type: "radar",
            data: {
                labels: ["GERP++ RS", "phyloP100way", "phastCons100way", "SiPhy", "LINSIGHT"],
                datasets: [{
                    data: [
                        asNumber(conservation["GERP++_RS"]) || 0,
                        asNumber(conservation.phyloP100way_vertebrate) || 0,
                        asNumber(conservation.phastCons100way_vertebrate) || 0,
                        asNumber(conservation.SiPhy_29way_logOdds) || 0,
                        asNumber(conservation.integrated_fitCons_score) || 0,
                    ],
                    backgroundColor: "rgba(15, 118, 110, 0.18)",
                    borderColor: "#0f766e",
                }],
            },
            options: { plugins: { legend: { display: false } } },
        });
    }

    async function batchSubmit(formData) {
        const progress = document.getElementById("batch-progress");
        let progressValue = 0;
        const timer = setInterval(() => {
            progressValue = Math.min(progressValue + 12, 85);
            if (progress) progress.style.width = `${progressValue}%`;
        }, 120);

        try {
            const data = await postForm("/api/v1/batch", formData);
            currentBatchPayload = data;
            localStorage.setItem("annotatemissenseBatch", JSON.stringify(data));
            if (progress) progress.style.width = "100%";
            window.location.href = "/batch_results";
            return data;
        } finally {
            clearInterval(timer);
        }
    }

    async function compareVariants(variantsArray) {
        const payload = await postJson("/api/v1/compare", { variants: variantsArray });

        const tableBody = document.getElementById("compareTableBody");
        const summary = document.getElementById("compareSummary");
        if (tableBody && summary) {
            tableBody.innerHTML = "";
            payload.variants.forEach((row) => {
                const tr = document.createElement("tr");
                tr.innerHTML = `
                    <td class="px-3 py-2">${row.chr}:${row.pos} ${row.ref}>${row.alt}</td>
                    <td class="px-3 py-2">${safeValue(row.gene)}</td>
                    <td class="px-3 py-2">${fmtScore(row.scores.annotatemissense_score)}</td>
                    <td class="px-3 py-2">${safeValue(row.scores.annotatemissense_pred)}</td>
                    <td class="px-3 py-2">${safeValue(row.scores.AlphaMissense)}</td>
                    <td class="px-3 py-2">${safeValue(row.scores.CADD_phred)}</td>
                    <td class="px-3 py-2">${safeValue(row.scores.REVEL_score)}</td>
                    <td class="px-3 py-2">${safeValue(row.scores.MetaRNN_score)}</td>
                `;
                tableBody.appendChild(tr);
            });
            summary.textContent = `${safeValue(payload.comparison.highest_annotatemissense_score)} has the highest AnnotateMissense score. Damaging votes per variant (${payload.comparison.tools}): ${payload.variants.map((r) => `${r.chr}:${r.pos} ${r.tool_agreement.damaging_votes}/${r.tool_agreement.tools_with_scores}`).join(" · ")}`;
        }

        const output = document.getElementById("compareOutput");
        if (output) output.textContent = JSON.stringify(payload, null, 2);
        return payload;
    }

    async function downloadCurrentVariantFullJson() {
        if (!currentSingleResult || !currentSingleResult.variant) {
            throw new Error("No variant loaded");
        }

        const v = currentSingleResult.variant;
        const params = new URLSearchParams({
            chrom: String(v.chr || ""),
            pos: String(v.pos || ""),
            ref: String(v.ref || ""),
            alt: String(v.alt || ""),
        });

        const response = await fetch(`/api/v1/variant/full?${params.toString()}`);
        const data = await response.json();
        if (!response.ok || data.status === "error") {
            throw new Error(data.message || "Unable to download full variant record");
        }

        const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = `variant_${String(v.chr).replace(/[^a-zA-Z0-9]/g, "_")}_${v.pos}_${v.ref}_${v.alt}_full.json`;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }

    function fillExample(chr, pos, ref, alt) {
        const chrom = document.getElementById("chrom");
        const p = document.getElementById("pos");
        const r = document.getElementById("ref");
        const a = document.getElementById("alt");
        if (!chrom || !p || !r || !a) return;
        chrom.value = chr;
        p.value = pos;
        r.value = ref;
        a.value = alt;
        singleVariantSearch(chr, pos, ref, alt).catch(() => {});
    }

    async function copyToClipboard(text) {
        if (!text) return;
        await navigator.clipboard.writeText(text);
    }

    function downloadResults(format) {
        const rows = currentBatchPayload?.results || currentGeneResults || [];
        if (!rows.length) return;
        const headers = ["chr", "pos", "ref", "alt", "status", "gene", "aa_change", "annotatemissense_score", "annotatemissense_confidence", "annotatemissense_pred", "AlphaMissense", "CADD_phred", "REVEL_score", "MetaRNN_score", "CLNSIG"];

        let content = "";
        let mime = "text/plain";
        let extension = format;

        if (format === "json") {
            content = JSON.stringify(rows, null, 2);
            mime = "application/json";
        } else {
            const sep = format === "tsv" ? "\t" : ",";
            const lines = [headers.join(sep)];
            rows.forEach((row) => {
                lines.push(headers.map((h) => {
                    const value = row[h] === null || row[h] === undefined ? "" : String(row[h]);
                    return sep === "," ? `"${value.replace(/"/g, '""')}"` : value;
                }).join(sep));
            });
            content = lines.join("\n");
            mime = "text/plain";
        }

        const blob = new Blob([content], { type: mime });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = `annotatemissense_results.${extension}`;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }

    function activateHomeTab(tabName) {
        document.querySelectorAll(".home-tab-button").forEach((btn) => {
            btn.classList.toggle("active", btn.dataset.homeTab === tabName);
        });
        document.querySelectorAll(".home-tab-panel").forEach((panel) => {
            panel.classList.toggle("hidden", panel.dataset.homePanel !== tabName);
        });
    }

    function initHomeTabs() {
        document.querySelectorAll(".home-tab-button").forEach((btn) => {
            btn.addEventListener("click", () => activateHomeTab(btn.dataset.homeTab));
        });
    }

    function initHomeSingleSearch() {
        const form = document.getElementById("single-variant-form");
        if (!form) return;
        form.addEventListener("submit", (event) => {
            event.preventDefault();
            const fd = new FormData(form);
            singleVariantSearch(fd.get("chrom"), fd.get("pos"), fd.get("ref"), fd.get("alt")).catch(() => {});
        });

        document.querySelectorAll("[data-example]").forEach((chip) => {
            chip.addEventListener("click", () => {
                const [chr, pos, ref, alt] = chip.dataset.example.split("|");
                fillExample(chr, pos, ref, alt);
            });
        });
    }

    function initBatchUpload() {
        const form = document.getElementById("batch-upload-form");
        if (!form) return;

        const dropZone = document.getElementById("dropZone");
        const fileInput = document.getElementById("variants_file");
        if (dropZone && fileInput) {
            ["dragenter", "dragover"].forEach((evt) => {
                dropZone.addEventListener(evt, (e) => {
                    e.preventDefault();
                    dropZone.classList.add("drag-over");
                });
            });
            ["dragleave", "drop"].forEach((evt) => {
                dropZone.addEventListener(evt, (e) => {
                    e.preventDefault();
                    dropZone.classList.remove("drag-over");
                });
            });
            dropZone.addEventListener("drop", (e) => {
                const files = e.dataTransfer?.files;
                if (files && files[0]) {
                    fileInput.files = files;
                }
            });
        }

        form.addEventListener("submit", async (event) => {
            event.preventDefault();
            const fd = new FormData(form);
            try {
                await batchSubmit(fd);
            } catch (error) {
                showInlineMessage(document.getElementById("single-result-preview"), error.message, "error");
            }
        });
    }

    function initGeneSearch() {
        const form = document.getElementById("gene-search-form");
        if (form) {
            form.addEventListener("submit", (event) => {
                event.preventDefault();
                const geneName = (new FormData(form).get("gene_name") || "").toString().trim();
                if (!geneName) return;
                window.location.href = `/gene/${encodeURIComponent(geneName)}`;
            });
        }

        const genePageForm = document.getElementById("genePageSearchForm");
        if (genePageForm) {
            genePageForm.addEventListener("submit", (event) => {
                event.preventDefault();
                const geneName = document.getElementById("genePageInput")?.value?.trim();
                if (!geneName) return;
                window.location.href = `/gene/${encodeURIComponent(geneName)}`;
            });
        }
    }

    function initResultsPage() {
        const payloadTag = document.getElementById("results-payload");
        const queryTag = document.getElementById("results-query");
        if (!payloadTag || !queryTag) return;

        const payload = JSON.parse(payloadTag.textContent || "{}");
        const query = JSON.parse(queryTag.textContent || "{}");

        if (payload && payload.status === "success") {
            renderResults(payload);
        } else if (!JSON.parse(document.getElementById("results-error")?.textContent || '""') && query.chr && query.pos && query.ref && query.alt) {
            singleVariantSearch(query.chr, query.pos, query.ref, query.alt).catch(() => {
                document.getElementById("resultsFallback")?.classList.remove("hidden");
            });
        }

        document.getElementById("download-json")?.addEventListener("click", async () => {
            try {
                await downloadCurrentVariantFullJson();
            } catch (error) {
                console.error(error);
            }
        });

        document.getElementById("share-link")?.addEventListener("click", async () => {
            await copyToClipboard(window.location.href);
        });

        document.querySelectorAll(".score-toggle").forEach((btn) => {
            btn.addEventListener("click", () => {
                document.querySelectorAll(".score-toggle").forEach((b) => b.classList.remove("active"));
                btn.classList.add("active");
            });
        });

        const drawer = document.getElementById("compareDrawer");
        document.getElementById("openCompareDrawer")?.addEventListener("click", () => drawer?.classList.remove("hidden"));
        document.getElementById("closeCompareDrawer")?.addEventListener("click", () => drawer?.classList.add("hidden"));
        document.getElementById("runCompareButton")?.addEventListener("click", () => {
            const text = document.getElementById("compareInput")?.value || "";
            const rows = text.split(/\r?\n/).map((v) => v.trim()).filter(Boolean);
            compareVariants(rows).catch((err) => {
                const output = document.getElementById("compareOutput");
                if (output) output.textContent = err.message;
            });
        });
    }

    function applyBatchFilterAndSort(results) {
        let rows = [...results];
        if (batchFilterState === "damaging") rows = rows.filter((r) => r.annotatemissense_pred === "D");
        if (batchFilterState === "clinvar") rows = rows.filter((r) => r.CLNSIG);
        if (batchFilterState === "missing") rows = rows.filter((r) => r.annotatemissense_score === null || r.annotatemissense_score === undefined);
        if (batchGeneFilter) rows = rows.filter((r) => String(r.gene || "").toLowerCase().includes(batchGeneFilter));

        rows.sort((a, b) => {
            const av = a[batchSortState.key];
            const bv = b[batchSortState.key];
            const an = asNumber(av);
            const bn = asNumber(bv);
            if (an !== null && bn !== null) return batchSortState.asc ? an - bn : bn - an;
            return batchSortState.asc
                ? String(av || "").localeCompare(String(bv || ""))
                : String(bv || "").localeCompare(String(av || ""));
        });
        return rows;
    }

    function renderBatchTable() {
        if (!currentBatchPayload) return;
        const body = document.getElementById("batchResultsBody");
        if (!body) return;

        const filtered = applyBatchFilterAndSort(currentBatchPayload.results || []);
        const pageSize = 20;
        const totalPages = Math.max(1, Math.ceil(filtered.length / pageSize));
        currentBatchPage = Math.min(currentBatchPage, totalPages);

        const start = (currentBatchPage - 1) * pageSize;
        const pageRows = filtered.slice(start, start + pageSize);

        body.innerHTML = "";
        pageRows.forEach((row) => {
            const tr = document.createElement("tr");
            tr.className = row.annotatemissense_pred === "D" ? "row-damaging" : row.annotatemissense_pred === "T" ? "row-tolerated" : "";
            tr.innerHTML = `
                <td class="px-3 py-2">${safeValue(row.variant || `${row.chr}:${row.pos} ${row.ref}>${row.alt}`)}</td>
                <td class="px-3 py-2">${safeValue(row.gene)}</td>
                <td class="px-3 py-2">${safeValue(row.annotatemissense_score)}</td>
                <td class="px-3 py-2">${safeValue(row.annotatemissense_pred)}</td>
                <td class="px-3 py-2">${safeValue(row.AlphaMissense)}</td>
                <td class="px-3 py-2">${safeValue(row.CADD_phred)}</td>
                <td class="px-3 py-2">${safeValue(row.REVEL_score)}</td>
                <td class="px-3 py-2">${safeValue(row.CLNSIG)}</td>
            `;
            body.appendChild(tr);

            const expand = document.createElement("tr");
            expand.className = "hidden";
            expand.innerHTML = `<td colspan="8" class="px-3 py-2 text-xs text-neutral">MetaRNN: ${safeValue(row.MetaRNN_score)} · gnomAD AF: ${safeValue(row.gnomad30_AF)}</td>`;
            tr.addEventListener("click", () => expand.classList.toggle("hidden"));
            body.appendChild(expand);
        });

        const pageInfo = document.getElementById("batchPageInfo");
        if (pageInfo) pageInfo.textContent = `Page ${currentBatchPage} of ${totalPages}`;
    }

    function renderBatchCharts() {
        if (!currentBatchPayload) return;
        const rows = currentBatchPayload.results || [];
        const damaging = rows.filter((r) => r.annotatemissense_pred === "D").length;
        const tolerated = rows.filter((r) => r.annotatemissense_pred === "T").length;

        drawChart("batchDonut", "batchDonutChart", {
            type: "doughnut",
            data: {
                labels: ["Damaging", "Tolerated", "Other"],
                datasets: [{ data: [damaging, tolerated, Math.max(0, rows.length - damaging - tolerated)], backgroundColor: ["#dc2626", "#16a34a", "#94a3b8"] }],
            },
        });

        const bins = [0, 0, 0, 0, 0];
        rows.forEach((r) => {
            const uq = asNumber(r.annotatemissense_score);
            if (uq === null) return;
            const idx = Math.min(4, Math.floor(uq * 5));
            bins[idx] += 1;
        });
        drawChart("batchHist", "batchHistogramChart", {
            type: "bar",
            data: {
                labels: ["0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"],
                datasets: [{ data: bins, backgroundColor: "#0f766e" }],
            },
            options: { plugins: { legend: { display: false } } },
        });

        const clinvarBuckets = { Pathogenic: 0, Benign: 0, VUS: 0, None: 0 };
        rows.forEach((r) => {
            const val = String(r.CLNSIG || "").toLowerCase();
            if (!val) clinvarBuckets.None += 1;
            else if (val.includes("pathogenic")) clinvarBuckets.Pathogenic += 1;
            else if (val.includes("benign")) clinvarBuckets.Benign += 1;
            else clinvarBuckets.VUS += 1;
        });
        drawChart("batchClinvar", "batchClinvarChart", {
            type: "bar",
            data: {
                labels: Object.keys(clinvarBuckets),
                datasets: [{ data: Object.values(clinvarBuckets), backgroundColor: ["#dc2626", "#16a34a", "#d97706", "#94a3b8"] }],
            },
            options: { plugins: { legend: { display: false } } },
        });
    }

    function initBatchResultsPage() {
        const table = document.getElementById("batchResultsTable");
        if (!table) return;
        const raw = localStorage.getItem("annotatemissenseBatch");
        if (!raw) return;
        currentBatchPayload = JSON.parse(raw);

        const summary = currentBatchPayload.summary || {};
        document.getElementById("sumSubmitted").textContent = safeValue(summary.submitted || currentBatchPayload.count);
        document.getElementById("sumFound").textContent = safeValue(summary.found);
        document.getElementById("sumNotFound").textContent = safeValue(summary.not_found);
        document.getElementById("sumDamaging").textContent = safeValue(summary.damaging);
        document.getElementById("sumTolerated").textContent = safeValue(summary.tolerated);

        table.querySelectorAll("th.sortable").forEach((th) => {
            th.addEventListener("click", () => {
                const key = th.dataset.key;
                if (batchSortState.key === key) batchSortState.asc = !batchSortState.asc;
                else batchSortState = { key, asc: true };
                renderBatchTable();
            });
        });

        document.querySelectorAll(".batch-filter").forEach((btn) => {
            btn.addEventListener("click", () => {
                batchFilterState = btn.dataset.filter;
                currentBatchPage = 1;
                renderBatchTable();
            });
        });

        document.getElementById("batchGeneSearch")?.addEventListener("input", (e) => {
            batchGeneFilter = e.target.value.trim().toLowerCase();
            currentBatchPage = 1;
            renderBatchTable();
        });

        document.getElementById("batchPrevPage")?.addEventListener("click", () => {
            currentBatchPage = Math.max(1, currentBatchPage - 1);
            renderBatchTable();
        });
        document.getElementById("batchNextPage")?.addEventListener("click", () => {
            currentBatchPage += 1;
            renderBatchTable();
        });

        document.querySelectorAll("[data-download-format]").forEach((btn) => {
            btn.addEventListener("click", () => downloadResults(btn.dataset.downloadFormat));
        });

        renderBatchCharts();
        renderBatchTable();
    }

    function initComparePage() {
        const compareSubmit = document.getElementById("compareSubmit");
        if (!compareSubmit) return;
        compareSubmit.addEventListener("click", () => {
            const text = document.getElementById("compareVariantsText")?.value || "";
            const rows = text.split(/\r?\n/).map((r) => r.trim()).filter(Boolean);
            compareVariants(rows).catch((error) => {
                const summary = document.getElementById("compareSummary");
                if (summary) summary.textContent = error.message;
            });
        });
    }

    function initApiDocs() {
        document.querySelectorAll(".copy-btn").forEach((btn) => {
            btn.addEventListener("click", async () => {
                const targetId = btn.dataset.copyTarget;
                const text = document.getElementById(targetId)?.innerText || "";
                await copyToClipboard(text);
                btn.textContent = "Copied";
                setTimeout(() => { btn.textContent = "Copy"; }, 1200);
            });
        });

        document.querySelectorAll("form.live-tester").forEach((form) => {
            form.addEventListener("submit", async (event) => {
                event.preventDefault();
                const block = form.closest(".api-endpoint-block");
                const output = block?.querySelector(".tester-output");
                const endpoint = block?.dataset.endpoint;
                const method = (block?.dataset.method || "GET").toUpperCase();
                if (!output || !endpoint) return;
                output.textContent = "Running...";
                try {
                    const response = await fetch(endpoint, { method, body: new FormData(form) });
                    const contentType = response.headers.get("content-type") || "";
                    if (contentType.includes("application/json")) {
                        const data = await response.json();
                        output.textContent = JSON.stringify(data, null, 2);
                    } else {
                        const text = await response.text();
                        output.textContent = `HTTP ${response.status}\n${text.slice(0, 2000)}`;
                    }
                } catch (error) {
                    output.textContent = error.message;
                }
            });
        });

        document.querySelectorAll(".run-get-tester").forEach((btn) => {
            btn.addEventListener("click", async () => {
                const block = btn.closest(".api-endpoint-block");
                const output = block?.querySelector(".tester-output");
                const endpoint = block?.dataset.endpoint;
                if (!output || !endpoint) return;
                output.textContent = "Running...";
                try {
                    const response = await fetch(endpoint);
                    const contentType = response.headers.get("content-type") || "";
                    if (contentType.includes("application/json")) {
                        const data = await response.json();
                        output.textContent = JSON.stringify(data, null, 2);
                    } else {
                        const text = await response.text();
                        output.textContent = `HTTP ${response.status}\n${text.slice(0, 2000)}`;
                    }
                } catch (error) {
                    output.textContent = error.message;
                }
            });
        });
    }

    function initGenePageData() {
        const dataTag = document.getElementById("gene-initial-json");
        if (!dataTag) return;
        currentGeneResults = JSON.parse(dataTag.textContent || "[]");

        document.querySelectorAll("[data-download-format]").forEach((btn) => {
            btn.addEventListener("click", () => downloadResults(btn.dataset.downloadFormat));
        });
    }

    function initAnnotatePage() {
        const form = document.getElementById("annotate-form");
        if (!form) return;
        const fileInput = document.getElementById("annotate_file");
        const textArea = document.getElementById("annotate_text");
        const message = document.getElementById("annotateMessage");
        const drop = document.getElementById("annotateDrop");

        if (drop && fileInput) {
            ["dragenter", "dragover"].forEach((evt) => drop.addEventListener(evt, (e) => { e.preventDefault(); drop.classList.add("drag-over"); }));
            ["dragleave", "drop"].forEach((evt) => drop.addEventListener(evt, (e) => { e.preventDefault(); drop.classList.remove("drag-over"); }));
            drop.addEventListener("drop", (e) => { if (e.dataTransfer?.files?.[0]) fileInput.files = e.dataTransfer.files; });
        }

        document.getElementById("loadExampleVcf")?.addEventListener("click", () => {
            textArea.value = (document.getElementById("example-vcf")?.textContent || "").trim() + "\n";
            if (fileInput) fileInput.value = "";
        });

        function buildForm() {
            const fd = new FormData();
            if (fileInput?.files?.[0]) fd.set("file", fileInput.files[0]);
            else if (textArea.value.trim()) fd.set("text", textArea.value);
            else throw new Error("Choose a VCF file or paste VCF text first");
            return fd;
        }

        function renderAnnotation(payload) {
            const summary = payload.summary || {};
            document.querySelectorAll("[data-sum]").forEach((el) => {
                const value = summary[el.dataset.sum];
                el.textContent = value === undefined ? "–" : Number(value).toLocaleString();
            });
            const warn = document.getElementById("annotateWarnings");
            if (warn) {
                warn.innerHTML = (payload.warnings || []).map((w) => `<div class="mb-2 rounded-lg border border-amber-200 bg-amber-50 p-2 text-xs text-amber-900">${escapeHtml(w)}</div>`).join("");
            }
            const body = document.getElementById("annotateBody");
            if (!body) return;
            const statusClass = { found: "text-success", no_score: "text-warning", not_in_release: "text-warning", not_snv: "text-neutral", invalid: "text-danger" };
            body.innerHTML = (payload.results || []).slice(0, 500).map((r) => {
                const link = r.query
                    ? `<a class="text-primary underline" href="/results?${new URLSearchParams({ chr: r.query.split(":")[0], pos: r.query.split(":")[1], ref: r.query.split(":")[2], alt: r.query.split(":")[3] })}">${escapeHtml(r.chrom)}:${escapeHtml(r.pos)} ${escapeHtml(r.ref)}&gt;${escapeHtml(r.alt)}</a>`
                    : `${escapeHtml(r.chrom)}:${escapeHtml(r.pos)} ${escapeHtml(r.ref)}&gt;${escapeHtml(r.alt)}`;
                const q = r.query ? r.query.split(":") : null;
                const changed = q && (q[1] !== String(r.pos) || q[2] !== String(r.ref).toUpperCase() || q[3] !== String(r.alt).toUpperCase());
                const rowClass = r.annotatemissense_pred === "D" ? "row-damaging" : r.annotatemissense_pred === "T" ? "row-tolerated" : "";
                return `<tr class="${rowClass}">
                    <td class="px-3 py-2 font-mono text-xs">${r.line}</td>
                    <td class="px-3 py-2 font-mono text-xs">${link}${changed ? `<br><span class="text-neutral">→ ${escapeHtml(r.query)}</span>` : ""}</td>
                    <td class="px-3 py-2 text-xs font-semibold ${statusClass[r.status] || ""}" title="${escapeHtml(r.message || "")}">${escapeHtml(r.status)}</td>
                    <td class="px-3 py-2 font-mono">${fmtScore(r.annotatemissense_score)}</td>
                    <td class="px-3 py-2">${safeValue(r.annotatemissense_pred)}</td>
                    <td class="px-3 py-2">${escapeHtml(r.gene || "--")}</td>
                    <td class="px-3 py-2 font-mono text-xs">${escapeHtml(String(r.aa_change_refgene || "--").split(",")[0])}</td>
                    <td class="px-3 py-2 text-xs">${escapeHtml(r.clinvar_clnsig || "--")}</td>
                </tr>`;
            }).join("") || '<tr><td colspan="8" class="px-3 py-6 text-center text-neutral">No records.</td></tr>';
        }

        form.addEventListener("submit", async (event) => {
            event.preventDefault();
            try {
                showInlineMessage(message, "Annotating…");
                const payload = await postForm("/api/v1/annotate?format=json", buildForm());
                renderAnnotation(payload);
                showInlineMessage(message, `Done in ${payload.seconds}s. Download the annotated VCF or TSV with the buttons above.`);
            } catch (error) {
                showInlineMessage(message, escapeHtml(error.message), "error");
            }
        });

        document.querySelectorAll("[data-annotate-download]").forEach((btn) => {
            btn.addEventListener("click", async () => {
                const fmt = btn.dataset.annotateDownload;
                try {
                    const response = await fetch(`/api/v1/annotate?format=${fmt}`, { method: "POST", body: buildForm() });
                    if (!response.ok) {
                        const err = await response.json().catch(() => ({}));
                        throw new Error(err.message || `HTTP ${response.status}`);
                    }
                    const blob = await response.blob();
                    const name = (response.headers.get("Content-Disposition") || "").split("filename=")[1] || `variants.annotatemissense.${fmt}`;
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement("a");
                    link.href = url;
                    link.download = name.replace(/"/g, "");
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(url);
                } catch (error) {
                    showInlineMessage(message, escapeHtml(error.message), "error");
                }
            });
        });
    }

    function initDefaultResultFetchFromUrl() {
        const onResultsPage = document.getElementById("results-query");
        if (onResultsPage) return;
        const params = new URLSearchParams(window.location.search);
        if (!["chr", "pos", "ref", "alt"].every((k) => params.get(k))) return;
        if (document.getElementById("single-variant-form")) {
            singleVariantSearch(params.get("chr"), params.get("pos"), params.get("ref"), params.get("alt")).catch(() => {});
        }
    }

    document.addEventListener("DOMContentLoaded", () => {
        initHomeTabs();
        initHomeSingleSearch();
        initBatchUpload();
        initGeneSearch();
        initResultsPage();
        initBatchResultsPage();
        initComparePage();
        initApiDocs();
        initGenePageData();
        initAnnotatePage();
        initDefaultResultFetchFromUrl();
    });

    window.singleVariantSearch = singleVariantSearch;
    window.renderResults = renderResults;
    window.batchSubmit = batchSubmit;
    window.compareVariants = compareVariants;
    window.fillExample = fillExample;
    window.copyToClipboard = copyToClipboard;
    window.downloadResults = downloadResults;
})();
