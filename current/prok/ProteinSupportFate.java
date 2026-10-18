package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.TreeMap;

import fileIO.ByteStreamWriter;
import prot.ProteinTailAssigner;
import shared.KillSwitch;
import structures.ByteBuilder;

/** Observation only: follows already-supported stops through native selection.
 * Coordinates are one-based genomic intervals; overlaps are not causal claims.
 * @author Sayu
 */
public final class ProteinSupportFate {

	ProteinSupportFate(ByteStreamWriter out_){
		require(out_!=null, "Fate output requires an owned writer");
		out=out_;
	}

	/** Call once per completed contig; caller owns writer header, ordering and close. */
	synchronized void write(ArrayList<Entry> entries, ArrayList<Orf> selected,
			int maxOverlapSame, int maxOverlapOpposite){
		try{
			final ByteBuilder rows=render(entries, selected, maxOverlapSame, maxOverlapOpposite);
			for(Entry e : entries){require(seen.add(e.key()), "Duplicate supported stop across contigs: "+e.key());}
			out.print(rows);
		}catch(RuntimeException | Error failure){
			KillSwitch.assertDie("Protein-support fate output failed: "+failure);
			throw failure;
		}
	}

	public static ByteBuilder render(ArrayList<Entry> entries, ArrayList<Orf> selected,
			int maxOverlapSame, int maxOverlapOpposite){
		require(entries!=null && selected!=null, "Fate rendering requires completed entry and path lists");
		final IdentityHashMap<Orf,Boolean> chosen=new IdentityHashMap<Orf,Boolean>();
		for(Orf o : selected){require(chosen.put(o, Boolean.TRUE)==null, "Duplicate selected object in final path");}
		final ArrayList<Entry> sorted=new ArrayList<Entry>(entries);
		Collections.sort(sorted, ORDER);
		final ByteBuilder rows=new ByteBuilder();
		String previous=null;
		for(Entry e : sorted){
			require(!e.key().equals(previous), "Duplicate supported stop: "+e.key()); previous=e.key();
			e.append(rows, chosen, selected, maxOverlapSame, maxOverlapOpposite);
		}
		return rows;
	}

	public static final class Entry {
		public Entry(Orf longest, ProteinSupport.StopEvidence evidence_, float stopScore_, int minLen_, float minStartScore_){
			require(longest!=null && longest.type==ProkObject.CDS && longest.flipped()==0,
				"Supported-stop hook requires strand-oriented longest CDS before candidate enumeration");
			require(evidence_!=null && evidence_.supported(), "Fate entries require observed model-supported stop evidence");
			rawName=longest.scafName; contig=contigToken(rawName); strand=longest.strand; scaflen=longest.scaflen;
			start=strand==0 ? longest.start+1 : scaflen-longest.stop;
			end=strand==0 ? longest.stop+1 : scaflen-longest.start;
			terminal=strand==0 ? end : start;
			require((strand==0 || strand==1) && start>0 && end>=start && end<=scaflen,
				"Longest-stop genomic coordinates must fit the contig");
			evidence=evidence_; stopScore=stopScore_; minLen=minLen_; minStartScore=minStartScore_;
		}

		public void admitted(Orf orf){
			require(orf!=null && orf.proteinEvidence!=null && orf.proteinEvidence.stop==evidence,
				"Admission hook must follow annotation with this exact stop's evidence");
			add(orf, orf.proteinEvidence.nativeScore);
		}

		private void add(Orf orf, float nativeScore){
			require(!sealed, "Admission after retained-list sealing");
			require(orf.type==ProkObject.CDS && orf.strand==strand && orf.flipped()==strand
				&& orf.scaflen==scaflen && orf.scafName.equals(rawName)
				&& (strand==0 ? orf.stop+1 : orf.start+1)==terminal,
				"Admitted candidate must use genomic coordinates and the registered supported stop");
			require(!byIdentity.containsKey(orf), "Candidate admitted twice");
			final Candidate c=new Candidate(orf, nativeScore);
			byIdentity.put(orf, c); candidates.add(c);
		}

		public void rejected(Orf orf, String firstPredicate){
			final Candidate c=byIdentity.get(orf);
			require(!sealed && c!=null && c.reason==null && firstPredicate!=null && !firstPredicate.isEmpty(),
				"Pruning must name one first predicate for an admitted, not-yet-rejected candidate");
			c.reason=firstPredicate;
		}

		public void retained(ArrayList<Orf> orfs){
			require(!sealed && orfs!=null, "Retained list must seal each stop exactly once, including empty lists");
			for(Orf orf : orfs){
				final Candidate c=byIdentity.get(orf);
				require(c!=null && c.reason==null && !c.retained, "Retained candidate must be admitted once and not rejected");
				c.retained=true;
			}
			for(Candidate c : candidates){require(c.retained || c.reason!=null, "Admitted candidate has neither retained nor rejected disposition");}
			sealed=true;
		}

		String key(){return escaped(contig)+":"+strand+":"+terminal;}

		private void append(ByteBuilder row, IdentityHashMap<Orf,Boolean> chosen, ArrayList<Orf> path, int same, int opposite){
			require(sealed && rejectedStartScore>=0 && rejectedLength>=0 && rejectedLimit>=0,
				"Final fate requires sealed accounting and nonnegative admission rejection counters");
			int retained=0, selected=0, rejected=0;
			Candidate best=null, bestRetained=null;
			final TreeMap<String,Integer> reasons=new TreeMap<String,Integer>();
			for(Candidate c : candidates){
				if(best==null || better(c, best)){best=c;}
				if(c.retained){retained++; if(bestRetained==null || better(c, bestRetained)){bestRetained=c;}}
				if(c.reason!=null){rejected++; final Integer n=reasons.get(c.reason); reasons.put(c.reason, n==null ? 1 : n+1);}
				if(chosen.containsKey(c.orf)){require(c.retained, "Final path selected a pruned candidate"); selected++;}
			}
			for(Orf o : path){
				if(o.type==ProkObject.CDS && o.scafName.equals(rawName) && o.strand==strand
					&& (strand==0 ? o.stop+1 : o.start+1)==terminal){
					require(byIdentity.containsKey(o), "Selected supported stop has an unobserved candidate identity");
				}
			}
			require(candidates.size()==retained+rejected && selected<=retained,
				"Stop accounting must conserve admitted=rejected+retained with selected a retained subset");
			final String fate=selected>0 ? "selected" : retained>0 ? "dp_not_selected" : candidates.isEmpty() ? "no_starts" : "soft_pruned";
			final ProteinTailAssigner.Hit h=evidence.hit;
			row.append("support-fate-v1").tab().append(key()).tab().append(escaped(contig)).tab().append(strand).tab().append(terminal)
				.tab().append(start).tab().append(end).tab().append(h.packIndex).tab().append(h.rank).tab().append(h.familyId)
				.tab().append(escaped(h.representative)).tab().append(h.identity, 8).tab().append(h.score64).tab().append(evidence.length)
				.tab().append(minLen).tab().appendSlow(minStartScore).tab().appendSlow(stopScore).tab().append(candidates.size())
				.tab().append(rejectedStartScore).tab().append(rejectedLength).tab().append(rejectedLimit)
				.tab().append(rejected).tab().append(retained).tab().append(selected).tab().append(fate).tab();
			if(reasons.isEmpty()){row.append('.');}
			else{boolean first=true; for(String reason : reasons.keySet()){
				if(!first){row.append(',');} first=false; row.append(escaped(reason)).append(':').append(reasons.get(reason).intValue());
			}}
			row.tab();
			if(best==null){row.append(".\t.\t.\t.\t.\t.\t.\t.");}
			else{row.append(best.start).tab().append(best.end).tab().appendSlow(best.rawStart).tab().appendSlow(best.rawStop)
				.tab().appendSlow(best.inner).tab().appendSlow(best.nativeScore).tab().appendSlow(best.bonus).tab().appendSlow(best.total);}
			final Candidate basis=bestRetained!=null ? bestRetained : best;
			row.tab().append(bestRetained!=null ? "best_retained" : best!=null ? "best_admitted" : "longest")
				.tab().append(basis==null ? start : basis.start).tab().append(basis==null ? end : basis.end)
				.tab().append(same).tab().append(opposite).tab();
			final TreeMap<String,String> overlaps=new TreeMap<String,String>();
			final int left=basis==null ? start : basis.start, right=basis==null ? end : basis.end;
			for(Orf o : path){
				if(!o.scafName.equals(rawName) || byIdentity.containsKey(o)){continue;}
				require(o.flipped()==o.strand, "Final-path overlap features must use genomic coordinates");
				final int overlap=Math.min(right, o.stop+1)-Math.max(left, o.start+1)+1;
				if(overlap<=0){continue;}
				final String id=escaped(contig)+":"+o.strand+":"+(o.start+1)+":"+(o.stop+1)+":"+o.type;
				final int limit=o.strand==strand ? same : opposite;
				require(overlaps.put(id, id+":"+overlap+":"+(overlap>limit ? 1 : 0))==null,
					"Duplicate selected feature key in overlap report: "+id);
			}
			if(overlaps.isEmpty()){row.append('.');}
			else{boolean first=true; for(String v : overlaps.values()){if(!first){row.append(',');} first=false; row.append(v);}}
			row.nl();
		}

		public int rejectedStartScore, rejectedLength, rejectedLimit;
		final String contig, rawName;
		final int strand, scaflen, start, end, terminal, minLen;
		final float stopScore, minStartScore;
		final ProteinSupport.StopEvidence evidence;
		private boolean sealed;
		private final ArrayList<Candidate> candidates=new ArrayList<Candidate>();
		private final IdentityHashMap<Orf,Candidate> byIdentity=new IdentityHashMap<Orf,Candidate>();
	}

	private static final class Candidate {
		Candidate(Orf o, float nativeScore_){
			assert(o!=null) : "Candidate snapshot follows validated admission";
			orf=o; start=o.start+1; end=o.stop+1; rawStart=o.startScore; rawStop=o.stopScore;
			inner=o.averageKmerScore(); nativeScore=nativeScore_; bonus=o.proteinBonus; total=o.orfScore;
		}
		final Orf orf;
		final int start, end;
		final float rawStart, rawStop, inner, nativeScore, bonus, total;
		String reason;
		boolean retained;
	}

	private static boolean better(Candidate a, Candidate b){
		assert(a!=null && b!=null) : "Best-candidate tie breaking needs two snapshots";
		return a.total>b.total || (a.total==b.total && (a.start<b.start || (a.start==b.start && a.end<b.end)));
	}
	private static String escaped(String s){
		require(s!=null, "Diagnostic identifiers cannot be null");
		return s.replace("%", "%25").replace(":", "%3A").replace(",", "%2C")
			.replace("\t", "%09").replace("\n", "%0A").replace("\r", "%0D");
	}
	/** Same first-space/tab identifier convention as Orf.appendGff; never rename caller objects. */
	public static String contigToken(String raw){
		require(raw!=null, "FASTA scaffold name cannot be null");
		int end=0;
		while(end<raw.length() && raw.charAt(end)!=' ' && raw.charAt(end)!='\t'){end++;}
		require(end>0, "FASTA scaffold name must have a nonempty emitted GFF identifier");
		return raw.substring(0, end);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalStateException(why);}}
	private static final Comparator<Entry> ORDER=new Comparator<Entry>(){
		@Override public int compare(Entry a, Entry b){
			int c=a.contig.compareTo(b.contig); if(c==0){c=Integer.compare(a.strand, b.strand);}
			return c==0 ? Integer.compare(a.terminal, b.terminal) : c;
		}
	};
	public static final String HEADER="schema\tstop_id\tcontig\tstrand\tterminal\tlongest_start\tlongest_end\tpack_index\tmodel_rank\tfamily_id\trepresentative\tidentity\tscore64\ttail_length\tmin_length\tmin_start_score\traw_stop_score\tadmitted\trejected_start_score\trejected_length\trejected_limit\tsoft_rejected\tretained\tselected\tfate\tfirst_pruning_predicates\tbest_start\tbest_end\tbest_raw_start\tbest_raw_stop\tbest_inner\tbest_native\tbest_bonus\tbest_total\toverlap_basis\toverlap_start\toverlap_end\tmax_overlap_same\tmax_overlap_opposite\tobserved_selected_overlaps";
	private final ByteStreamWriter out;
	private final HashSet<String> seen=new HashSet<String>();
}
