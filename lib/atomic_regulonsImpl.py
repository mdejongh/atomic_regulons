#BEGIN_HEADER
#END_HEADER


class atomic_regulons:
    '''
    Module Name:
    atomic_regulons

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
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
            raise ValueError('Method compute_atomic_regulons return value ' +
                             'atomic_regulons is not type list as required.')
        if not isinstance(feature_calls, list):
            raise ValueError('Method compute_atomic_regulons return value ' +
                             'feature_calls is not type list as required.')
        if not isinstance(ar_calls, list):
            raise ValueError('Method compute_atomic_regulons return value ' +
                             'ar_calls is not type list as required.')
        # return the results
        return [atomic_regulons, feature_calls, ar_calls]

    def compute_atomic_regulons_CDS(self, genome_id):
        # self.ctx is set by the wsgi application class
        # return variables are: atomic_regulons, feature_calls, ar_calls
        #BEGIN compute_atomic_regulons_CDS
        #END compute_atomic_regulons_CDS

        #At some point might do deeper type checking...
        if not isinstance(atomic_regulons, list):
            raise ValueError('Method compute_atomic_regulons_CDS return value ' +
                             'atomic_regulons is not type list as required.')
        if not isinstance(feature_calls, list):
            raise ValueError('Method compute_atomic_regulons_CDS return value ' +
                             'feature_calls is not type list as required.')
        if not isinstance(ar_calls, list):
            raise ValueError('Method compute_atomic_regulons_CDS return value ' +
                             'ar_calls is not type list as required.')
        # return the results
        return [atomic_regulons, feature_calls, ar_calls]

    def compute_atomic_regulons_expressionServices(self, genome_id, sample_ids):
        # self.ctx is set by the wsgi application class
        # return variables are: atomic_regulons, feature_calls, ar_calls
        #BEGIN compute_atomic_regulons_expressionServices
        #END compute_atomic_regulons_expressionServices

        #At some point might do deeper type checking...
        if not isinstance(atomic_regulons, list):
            raise ValueError('Method compute_atomic_regulons_expressionServices return value ' +
                             'atomic_regulons is not type list as required.')
        if not isinstance(feature_calls, list):
            raise ValueError('Method compute_atomic_regulons_expressionServices return value ' +
                             'feature_calls is not type list as required.')
        if not isinstance(ar_calls, list):
            raise ValueError('Method compute_atomic_regulons_expressionServices return value ' +
                             'ar_calls is not type list as required.')
        # return the results
        return [atomic_regulons, feature_calls, ar_calls]
