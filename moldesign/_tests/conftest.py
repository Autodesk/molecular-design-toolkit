

def pytest_itemcollected(item):
    marks = getattr(item.module, '__PYTEST_MARK__', None)
    if marks is None:
        return
    if isinstance(marks, str):
        marks = [marks]
    for mark in marks:
        item.add_marker(mark)


# TODO: nicer output strings for git commit status
# see https://docs.pytest.org/en/latest/example/simple.html#post-process-test-reports-failures
#@pytest.hookimpl(tryfirst=True, hookwrapper=True)
#def pytest_runtest_makereport(item, call):i
#    pass
# Also possibly useful:  item.add_report_section

