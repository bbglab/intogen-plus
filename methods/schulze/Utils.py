def join_dictionaries(list_dictionaries):
    '''
    Given a list of dictionaries join them and ouputs an unique one
    :param list_dictionaries: list of dictionaries
    :return: joint dictionary
    '''
    result = {}
    for d in list_dictionaries:

        for cancer in d.keys():
            if not cancer in result:
                result[cancer] = {}
            for method in d[cancer]:

                result[cancer][method] = dict(d[cancer][method])

    return result
