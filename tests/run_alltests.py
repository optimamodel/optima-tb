tests = ['project', 'model', 'programs']
prefix = 'test_'
suffix = '.py'

for test in tests:
    test_type  = prefix + test + suffix
    execfile(test_type)