# !/bin/bash

echo "Uncompressing the assignment.."
tar -zxvf assignment.tar.gz
cd assignment/

echo "Installing dependencies.."
pip install -r requirements.txt
echo "Setup complete.."
