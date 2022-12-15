rm -r dist ;
python setup.py sdist bdist_wheel ;
pip install dist/FACS_Samplin-0.1.0-py3-none-any.whl ; #--force-reinstall;
#if twine check dist/* ; then
#  if [ "$1" = "--test" ] ; then
#    twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#  else
#    twine upload dist/* ;
#  fi
#fi