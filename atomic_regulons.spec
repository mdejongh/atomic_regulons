module atomic_regulons {
    /* genome id needs to be KBase interpretable */
    typedef string genome_id;

    /* table of expression levels for features in samples */
    typedef structure {
	list<string> sample_names;
	mapping<string feature_id, list<string> expression_levels> expression_vectors;
    } ExpressionValues;

    /* atomic regulon has an id and a set of features */
    typedef structure {
	string ar_id;
	list<string> feature_ids;
    } AtomicRegulon;

    /* on/off/unknown call for a feature in one sample */
    typedef structure {
	string sample_name;
	string feature_id;
	int on_off_unknown;
    } FeatureOnOffCall;

    /* on/off/unknown call for an AtomicRegulon in one sample */
    typedef structure {
	string sample_name;
	string ar_id;
	int on_off_unknown;
    } AtomicRegulonOnOffCall;
    

    /* input a list of expression values for a genome, compute atomic regulons */
    funcdef compute_atomic_regulons(ExpressionValues expression_values, string genome_id) returns (list<AtomicRegulon> atomic_regulons, list<FeatureOnOffCall> feature_calls, list<AtomicRegulonOnOffCall> ar_calls);
};