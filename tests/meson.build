test_names = [
    'centroid',
    'ellipse',
#    'gaussian',
]
foreach test_name : test_names
    test = executable(
        'gauss2_test_' + test_name,
        'test_' + test_name + '.cc',
        include_directories : inc,
        link_with : gauss2dlib,
    )
    test(test_name, test)
endforeach