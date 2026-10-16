package prok;

import java.util.HashMap;
import java.nio.charset.StandardCharsets;

import fileIO.ByteFile;
import parse.LineParser1;
import prot.ProteinFamilyAssigner;

/** Optional, roster-bound selected labels; public search tools are build-time only.
 * @author Keqing
 */
final class FamilyLabels {

	static FamilyLabels load(final String path, final ProteinFamilyAssigner library){
		final int n=library.familyCount(); final int[] ids=new int[n]; final String[] reps=new String[n];
		for(int i=0; i<n; i++){ids[i]=library.familyId(i); reps[i]=library.representative(i);}
		return load(path, ids, reps, library.identityPin(), library.consensusPin());
	}

	/** Small metadata-only entry point also permits tests without loading a model. */
	static FamilyLabels load(final String path, final int[] ids, final String[] representatives,
			final String identityPin, final String consensusPin){
		require(ids!=null && representatives!=null && ids.length==representatives.length, "Label reference dimensions differ");
		final HashMap<String,String> meta=new HashMap<String,String>();
		final HashMap<String,Integer> columns=new HashMap<String,Integer>();
		final Label[] labels=new Label[ids.length];
		final ByteFile input=ByteFile.makeByteFile(path, false); final LineParser1 lp=new LineParser1('\t');
		String version=null;
		int schema=0;
		int width=0, row=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				row++; if(line.length==0){continue;} lp.set(line);
				if(line[0]=='#'){
					require(width==0 && lp.terms()==2, "Label metadata must precede the TSV header and have two fields: row "+row);
					final String key=lp.parseString(0).substring(1), value=lp.parseString(1);
					require(!value.isEmpty() && meta.putIfAbsent(key, value)==null, "Empty or duplicate label metadata: "+key);
					continue;
				}
				if(width==0){
					schema=Integer.parseInt(meta.get("schema_version"));
					require((schema==1 || schema==2) && "bbgap_family_labels".equals(meta.get("kind")), "Unsupported label schema or kind");
					require(identityPin.equals(meta.get("family_identity_sha80")), "Label identity pin differs from loaded family roster");
					require(consensusPin.equals(meta.get("consensus_sha80")), "Label consensus pin differs from loaded sequences");
					version=meta.get("label_version"); require(version!=null && version.matches("[A-Za-z0-9._-]+"), "Label version must be a nonempty ASCII release token");
					width=lp.terms();
					for(int i=0; i<width; i++){
						final String name=lp.parseString(i);
						require(!name.isEmpty() && columns.putIfAbsent(name, i)==null, "Empty or duplicate label column: "+name);
					}
					for(String key : REQUIRED){require(columns.containsKey(key), "Missing label column: "+key);}
					if(schema==2){
						for(String key : REQUIRED_V2){require(columns.containsKey(key), "Missing label column: "+key);}
					}
					continue;
				}
				require(lp.terms()==width, "Wrong label column count at row "+row);
				final int rank=lp.parseInt(columns.get("active_index"));
				require(rank>=0 && rank<ids.length && labels[rank]==null, "Duplicate or unknown label rank at row "+row);
				require(lp.parseInt(columns.get("family_id"))==ids[rank]
					&& lp.termEquals(representatives[rank], columns.get("rep_id")), "Label family ID/representative mismatch at row "+row);
				require(lp.termEquals(version, columns.get("label_version")), "Mixed label releases at row "+row);
				final String decision=lp.parseString(columns.get("decision_class"));
				final String product=parseUtf8(lp, columns.get("product"));
				final boolean provisional=isProvisional(decision);
				final boolean named=decision.equals("specific_full_length") || decision.equals("swissprot_full_length")
					|| decision.equals("domain_supported_broad_name") || decision.equals("homology_supported_broad_name")
					|| provisional;
				require(named || decision.equals("conflicting_evidence") || decision.equals("partial_match") || decision.equals("unnamed"),
					"Unknown label decision class at row "+row+": "+decision);
				require(!named || !product.isEmpty(), "Named label has no product at row "+row);
				require(!provisional || schema==2, "Provisional labels require schema 2 at row "+row);
				final String display=schema==2 ? parseUtf8(lp, columns.get("display_label")) : null;
				final String source=schema==2 ? parseUtf8(lp, columns.get("source_methods")) : null;
				final String support=schema==2 ? parseUtf8(lp, columns.get("support_tag")) : null;
				final String coverage=schema==2 ? parseUtf8(lp, columns.get("coverage_tag")) : null;
				final String score=schema==2 ? parseUtf8(lp, columns.get("score_tag")) : null;
				final String productText=schema==2 ? parseUtf8(lp, columns.get("product_text_tag")) : null;
				final String review=schema==2 ? parseUtf8(lp, columns.get("review_tag")) : null;
				final String hit=schema==2 ? parseUtf8(lp, columns.get("selected_hit_id")) : null;
				require(!provisional || (!display.isEmpty() && !source.isEmpty()
					&& !support.isEmpty() && !coverage.isEmpty() && !score.isEmpty()
					&& !productText.isEmpty() && !review.isEmpty() && !hit.isEmpty()),
					"Provisional schema-2 labels require complete evidence fields at row "+row);
				labels[rank]=new Label(decision, named ? product : null,
					display, source, support, coverage, score, productText, review, hit);
			}
		}finally{require(!input.close(), "Label input failed: "+path);}
		require(width>0, "Missing label TSV header");
		return new FamilyLabels(version, labels);
	}

	private FamilyLabels(String version_, Label[] labels_){version=version_; labels=labels_;}
	Label get(int rank){return labels[rank];}
	private static String parseUtf8(final LineParser1 lp, final int column){
		return new String(lp.parseByteArray(column), StandardCharsets.UTF_8);
	}
	private static boolean isProvisional(final String decision){
		return decision.equals("provisional_full") || decision.equals("weak_provisional")
			|| decision.equals("weak_provisional_placeholder");
	}
	static final class Label {
		Label(String decision_, String product_, String display_, String source_, String support_,
				String coverage_, String score_, String productText_, String review_,
				String hit_){
			decision=decision_;
			product=product_;
			display=display_;
			source=source_;
			support=support_;
			coverage=coverage_;
			score=score_;
			productText=productText_;
			review=review_;
			hit=hit_;
		}
		final String decision, product;
		final String display, source, support, coverage, score, productText, review, hit;
	}
	private static void require(boolean ok, String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	final String version;
	private final Label[] labels;
	private static final String[] REQUIRED={"family_id", "active_index", "rep_id", "label_version", "decision_class", "product"};
	private static final String[] REQUIRED_V2={"display_label", "source_methods", "support_tag", "coverage_tag", "score_tag", "product_text_tag", "review_tag", "selected_hit_id"};
}
