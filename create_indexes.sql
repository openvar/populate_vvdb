ALTER TABLE stableGeneIds ADD INDEX hgnc_id_index (hgnc_id);
ALTER TABLE stableGeneIds ADD INDEX hgnc_symbol_index (hgnc_symbol);
ALTER TABLE transcript_info ADD INDEX refSeqID_index (refSeqID);
ALTER TABLE transcript_info ADD INDEX utaSymbol_index (utaSymbol);
ALTER TABLE transcript_info ADD INDEX hgncSymbol_index_index (hgncSymbol);
ALTER TABLE LRG_proteins ADD INDEX RefSeqProteinID_index (RefSeqProteinID);
ALTER TABLE LRG_transcripts ADD INDEX RefSeqTranscriptID_index (RefSeqTranscriptID);
