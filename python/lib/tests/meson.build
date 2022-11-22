test_names = [
    'pyimage', # commented out until it can pass locally
]
pymod = import('python')
python = pymod.find_installation(
    'python3',
    required : true,
)
foreach test_name : test_names
    test = executable(
        'gauss2_python_test_' + test_name,
        'test_' + test_name + '.cc',
        dependencies: [ # List of dependencies
            python.dependency(embed: true), # Add Python.h as dependency
            gauss2d_dep,
            pybind11_dep,
        ],
        include_directories : inc,
    )
    test(test_name, test)
endforeach
