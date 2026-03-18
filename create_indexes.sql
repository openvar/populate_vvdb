ALTER TABLE stableGeneIds ADD INDEX hgnc_id_index (hgnc_id);
ALTER TABLE refSeqGene_loci ADD INDEX refSeqGeneID_index (refSeqGeneID);
ALTER TABLE transcript_info ADD INDEX refSeqID_index (refSeqID);
ALTER TABLE transcript_info ADD INDEX utaSymbol_index (utaSymbol);
ALTER TABLE transcript_info ADD INDEX hgncSymbol_index (hgncSymbol);


ALTER TABLE stableGeneIds DROP INDEX hgnc_id_index;
ALTER TABLE refSeqGene_loci DROP INDEX refSeqGeneID_index;
ALTER TABLE transcript_info DROP INDEX refSeqID_index;
ALTER TABLE transcript_info DROP INDEX utaSymbol_index;
ALTER TABLE transcript_info DROP INDEX hgncSymbol_index;