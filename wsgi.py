"""WSGI entry point: gunicorn wsgi:app  /  waitress-serve wsgi:app"""
from app import create_app

app = create_app()
