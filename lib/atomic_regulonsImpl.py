#BEGIN_HEADER
#END_HEADER

'''

Module Name:
atomic_regulons

Module Description:


'''
class atomic_regulons:

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    def __init__(self, config): #config contains contents of config file in hash or 
                                #None if it couldn't be found
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def compute_atomic_regulons(self, expression_values, genome_id):
        # self.ctx is set by the wsgi application class
        # return variables are: atomic_regulons, feature_calls, ar_calls
        #BEGIN compute_atomic_regulons
        #END compute_atomic_regulons

        #At some point might do deeper type checking...
        if not isinstance(atomic_regulons, list):
            raise ValueError('Method compute_atomic_regulons return value atomic_regulons is not type list as required.')
        if not isinstance(feature_calls, list):
            raise ValueError('Method compute_atomic_regulons return value feature_calls is not type list as required.')
        if not isinstance(ar_calls, list):
            raise ValueError('Method compute_atomic_regulons return value ar_calls is not type list as required.')
        # return the results
        return [ atomic_regulons, feature_calls, ar_calls ]
        
