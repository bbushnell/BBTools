package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Reports independent CM evidence for prediction-derived review contexts.
 * Every raw hit survives. No annotation edit, automatic admission or FP verdict.
 * @author Brian, Raiden
 */
public final class NcrnaReviewCmEvidence {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaReviewCmEvidence.class, false);
		String loci=null, hits=null, out=null, model=null, family=null;int expected=-1, clen=-1;double ga=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("loci")){loci=value;}else if(key.equals("hits")){hits=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("modelid")){model=value;}else if(key.equals("family")){family=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else if(key.equals("clen")){clen=Integer.parseInt(value);}
			else if(key.equals("gabits")){ga=Double.parseDouble(value);}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(loci!=null && hits!=null && out!=null && expected>0, "Bind prediction metadata, CM table and union denominator");
		final RrnaCmTruthBuilder.Criteria criteria=new RrnaCmTruthBuilder.Criteria(model, family, clen, clen, ga);
		Shared.setThreads(1);
		try{new NcrnaReviewCmEvidence().run(loci, hits, out, expected, criteria);}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(final String metadata, final String table, final String out, final int expected, final RrnaCmTruthBuilder.Criteria criteria){
		final ArrayList<Locus> loci=new ArrayList<Locus>();final ObjectMap<String,Locus> ids=new ObjectMap<String,Locus>(String.class, Locus.class);
		final ByteFile in=RrnaResourceIO.open(metadata);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaReviewUnion.LOCUS_HEADER, metadata);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==14 && p.termEquals(criteria.family, 5) && p.termEquals("PREDICTION_UNRESOLVED", 13), "Only the declared prediction-derived union may be evaluated");
				final Locus l=new Locus(p);
				require(!ids.contains(l.id), "Repeated review ID would merge independent evidence: "+l.id);ids.put(l.id, l);loci.add(l);
			}
		}finally{require(!in.close(), "Review metadata read failed");}
		require(loci.size()==expected, "Review union denominator changed");
		final ArrayList<RrnaCmTruthBuilder.Hit> hits=RrnaCmTruthBuilder.readHits(new RrnaCmTruthBuilder.Genome("review", "contexts", table, null), criteria);
		for(RrnaCmTruthBuilder.Hit h:hits){
			final Locus l=ids.get(h.seqid);require(l!=null && Math.max(h.seqFrom, h.seqTo)<=l.length, "Every CM hit must belong inside one frozen review context");l.hits.add(h);
		}
		RrnaResourceIO.fresh(out, ".hits.tsv", ".loci.tsv", ".summary.tsv");
		final ByteStreamWriter raw=RrnaResourceIO.writer(out+".hits.tsv"), rows=RrnaResourceIO.writer(out+".loci.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int none=0, one=0, multiple=0, events=0;
		try{
			raw.println("review_id\tcm_line\tseqid\tgenomic_start1\tgenomic_stop1\tgenomic_strand\tmodel_from\tmodel_to\ttrunc\tbits\tevalue\tincluded\tsense_overlap_bp\tfull_ga_reciprocal_half");
			rows.println("review_id\traw_hits\tsense_overlapping_hits\tfull_ga_reciprocal_half_hits\tstatus");
			for(Locus l:loci){
				int overlapping=0, qualified=0;
				for(RrnaCmTruthBuilder.Hit h:l.hits){
					final int from=Math.min(h.seqFrom, h.seqTo), to=Math.max(h.seqFrom, h.seqTo), callLength=l.stop-l.start+1;
					final int overlap=h.strand=='+'?Math.max(0, Math.min(to, l.left+callLength)-Math.max(from, l.left+1)+1):0;
					final boolean pass=overlap>0 && 2L*overlap>=callLength && 2L*overlap>=to-from+1
						&& h.included && Double.parseDouble(h.score)>=criteria.ga && h.mdlFrom==1 && h.mdlTo==criteria.clen && h.trunc.equals("no");
					if(overlap>0){overlapping++;}if(pass){qualified++;}events++;
					final int a=l.strand=='+'?l.start-l.left+from-1:l.stop+l.left-to+1;
					final int z=l.strand=='+'?l.start-l.left+to-1:l.stop+l.left-from+1;
					final char strand=h.strand=='+'?l.strand:l.strand=='+'?'-':'+';
					require(a>0 && z>=a && z<=l.contigLength && z-a==to-from, "CM projection must preserve span within its original source contig");
					raw.print(b.clear().append(l.id).tab().append(h.line).tab().append(l.seqid).tab().append(a).tab().append(z).tab().append(strand)
						.tab().append(h.mdlFrom).tab().append(h.mdlTo).tab().append(h.trunc).tab().append(h.score).tab().append(h.evalue)
						.tab().append(h.included).tab().append(overlap).tab().append(pass).nl());
				}
				final String status=qualified==0?"NO_FULL_GA_OVERLAPPING_HIT":qualified==1?"ONE_FULL_GA_OVERLAPPING_HIT":"MULTIPLE_FULL_GA_OVERLAPPING_HITS";
				if(qualified==0){none++;}else if(qualified==1){one++;}else{multiple++;}
				rows.print(b.clear().append(l.id).tab().append(l.hits.size()).tab().append(overlapping).tab().append(qualified).tab().append(status).nl());
			}
			require(events==hits.size() && none+one+multiple==expected, "Independent evidence must conserve every raw hit and every predicted review locus");
			summary.println("loci\traw_hits\tno_qualifying_hit\tone_qualifying_hit\tmultiple_qualifying_hits\tinterpretation");
			summary.print(b.clear().append(expected).tab().append(events).tab().append(none).tab().append(one).tab().append(multiple)
				.append("\tCM_EVIDENCE_NOT_TRUTH_EDITS_OR_FP_VERDICTS\n"));
		}finally{RrnaResourceIO.finish(raw, rows, summary);}
		System.err.println("NCRNA_REVIEW_CM_EVIDENCE_PASS loci="+expected+" rawHits="+events+" oneFull="+one+" multipleFull="+multiple);
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Locus {
		Locus(LineParser1 p){
			id=p.parseString(0);seqid=p.parseString(4);final String s=p.parseString(6);
			require(s.equals("+") || s.equals("-"), "Prediction strand must be explicit");strand=s.charAt(0);
			start=p.parseInt(7);stop=p.parseInt(8);contigLength=p.parseInt(9);left=p.parseInt(10);final int right=p.parseInt(11);length=p.parseInt(12);
			require(start>0 && stop>=start && stop<=contigLength && left>=0 && right>=0 && (long)length==(long)stop-start+1+left+right,
				"Review metadata must preserve the original predicted span and actual sense flanks");
			require(strand=='+'?start-left>=1 && (long)stop+right<=contigLength:start-right>=1 && (long)stop+left<=contigLength,
				"Clipped context must remain inside the original contig on its declared strand");
		}
		final String id, seqid;final char strand;final int start, stop, contigLength, left, length;
		final ArrayList<RrnaCmTruthBuilder.Hit> hits=new ArrayList<RrnaCmTruthBuilder.Hit>();
	}
}
