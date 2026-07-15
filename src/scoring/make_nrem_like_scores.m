function nremLike = make_nrem_like_scores(analysisCode, cfg)
%MAKE_NREM_LIKE_SCORES Merge microarousals into NREM for bout extraction only.

nremLike = int16(analysisCode(:));
nremLike(nremLike == cfg.code_microarousal) = int16(cfg.code_nrem);
end
