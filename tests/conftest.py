import pytest
from app.app import app  # Make sure app is correctly imported from your application

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client
