@echo off

echo Uncompressing the assignment...
tar -zxvf assignment.tar.gz
cd assignment\

echo Installing dependencies...
cd code_index_query\code\
pip install -r requirements.txt

echo Setup complete...
pause
