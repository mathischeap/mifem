

"""Encrypt_password"""


# %% generate key
from cryptography.fernet import Fernet
key = Fernet.generate_key()
print(key)

# %% encrypt a password
key = b'3H148a_5zmo9b1FQCWTnte4NCUdwP8NRbgAq0eU6D7A='
cipher_suite = Fernet(key)
ciphered_text = cipher_suite.encrypt(b"ASDAS")   #required to be bytes
print(ciphered_text)

# %% recover a password
key = b'3H148a_5zmo9b1FQCWTnte4NCUdwP8NRbgAq0eU6D7A='

cipher_suite = Fernet(key)
ciphered_text = b'gAAAAABdRYkBbflLkADtKxwODgRXSDW5tfMcv7MhBxNb6WGmif3q6LWrDgsfN3TVz51Ybn3-VJAYfk2l7c0A9JYbH6NE9z5DbQ=='
unciphered_text = cipher_suite.decrypt(ciphered_text)
print(unciphered_text)