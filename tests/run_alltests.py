import unittest

testmodules = []
tests = ['project', 'programs']
prefix = 'test_'
suffix = '.py'

for test in tests: testmodules.append(prefix+test)

suite = unittest.TestSuite()

for t in testmodules:
    suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))
test_runner = unittest.TextTestRunner().run(suite)