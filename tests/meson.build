test_names = [
    'centroid',
    'coordinatesystem',
    'ellipse',
    'evaluate',
    'gaussian',
    'image',
]
foreach test_name : test_names
    test = executable(
        'gauss2d_test_' + test_name,
        'test_' + test_name + '.cc',
        include_directories : inc,
        link_with : gauss2dlib,
    )
    test(test_name, test)
endforeach
