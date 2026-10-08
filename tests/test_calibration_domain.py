from microc2compare.science import calibration_design_domain,calibration_key,design_matrix


def test_calibration_key_does_not_transfer_from_three_to_two_units(config):
    comparison={'comparison_id':'x','reference_condition':'C','case_condition':'T','design':'~condition','inference':'formal'}
    def samples(n):return [{'sample_id':c+str(i),'condition':c,'subject_id':c+str(i)} for c in ('C','T') for i in range(n)]
    three=calibration_design_domain(config,comparison,samples(3))
    comparison=dict(comparison,inference='preliminary')
    two=calibration_design_domain(config,comparison,samples(2))
    assert three['residual_df']==4 and two['residual_df']==2
    assert calibration_key(config,'loops',comparison,{},three)!=calibration_key(config,'loops',comparison,{},two)


def test_numeric_subject_labels_are_factors(config):
    comparison={'comparison_id':'x','reference_condition':'C','case_condition':'T','design':'~subject_id+condition','inference':'formal'}
    samples=[{'sample_id':c+str(i),'condition':c,'subject_id':str(i)} for c in ('C','T') for i in (1,2,3)]
    design=design_matrix(samples,comparison,config)
    assert len(design['columns'])==4 and design['residual_df']==2
