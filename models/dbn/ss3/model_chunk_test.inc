variable: change_state {
	type: discrete hidden cardinality BINARY_CARD;
	switchingparents: nil;
	conditionalparents: state_count_down(0) using DenseCPT("Fi_given_Di_CPT");
}

variable: state {
	type: discrete hidden cardinality STATE_CARD;
	switchingparents: change_state(-1) using mapping("1D_map");
	conditionalparents:
		state(-1) using DenseCPT("state_extension_CPT")
		| state(-1) using DenseCPT("state_transition_CPT");
}

variable: state_count_down {
	type: discrete hidden cardinality MAX_STATE_LENGTH_CARD;
	switchingparents: change_state(-1) using mapping("1D_map");
	conditionalparents:
		state_count_down(-1), state(0) using DenseCPT("exponential_tail_CPT")
		| state(0) using DenseCPT("length_CPT");
}

variable: state_class_history {
	type: discrete hidden cardinality STATE_CLASS_HISTORY_CARD;
	switchingparents: nil;
	conditionalparents: state_class_history(-1), state(0) using DeterministicCPT("update_history_CPT");
}

variable: aa_profile {
	type: continuous observed 0:19;
	weight: scale WEIGHT_AA_PROFILE;
	switchingparents: nil;
	conditionalparents: state_class_history(0) using mixture collection("ECcollection") mapping("tuple_map");
}

