ALTER TABLE stableGeneIds ADD INDEX hgnc_id_index (hgnc_id);
ALTER TABLE transcript_info ADD INDEX refSeqID_index (refSeqID);
ALTER TABLE transcript_info ADD INDEX utaSymbol_index (utaSymbol);
ALTER TABLE transcript_info ADD INDEX hgncSymbol_index_index (hgncSymbol);