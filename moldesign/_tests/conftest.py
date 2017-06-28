

def pytest_itemcollected(item):
    if hasattr(item.module, '__PYTEST_MARK__'):
        item.add_marker(item.module.__PYTEST_MARK__)


# TODO: nicer output strings for git commit status
# see https://docs.pytest.org/en/latest/example/simple.html#post-process-test-reports-failures
#@pytest.hookimpl(tryfirst=True, hookwrapper=True)
#def pytest_runtest_makereport(item, call):i
#    pass
# Also possibly useful:  item.add_report_section

