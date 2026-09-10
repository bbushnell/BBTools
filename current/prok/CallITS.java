package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import gff.GffLine;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Emits fungal ITS regions from a caller GFF containing 18S, 5.8S and LSU/28S rRNA rows
 * (the rrnafallback three-flank output: legacy 18S via 18s=t plus the r58/lsu scavenger
 * families). Implements the reviewed ITS emission contract
 * (rrnafallback plans/fungal_its_and_resource_handoff_20260909.md; acceptance authority =
 * Qiqi's fungal_its_synthetic_acceptance_fixtures_20260909.md and
 * its_interval_emitter_edge_case_review_20260909.md):
 *
 * <ul>
 * <li>PRIMARY output (Brian, 2026-09-09: "for most uses you want ITS1+5.8S+ITS2 as a full
 *   contiguous sequence") = fullITS bounded by the OUTER flanks: plus
 *   [18S.end+1, LSU.start-1], minus [LSU.end+1, 18S.start-1]. A 5.8S pivot is NOT required
 *   when the outer pairing is unambiguous; it assists pairing but is never mandatory.</li>
 * <li>Individual ITS1/ITS2 are secondary outputs, each conditional on a uniquely paired
 *   5.8S pivot: plus ITS1=[18S.end+1, 5.8S.start-1], ITS2=[5.8S.end+1, LSU.start-1];
 *   minus ITS1=[5.8S.end+1, 18S.start-1], ITS2=[LSU.end+1, 5.8S.start-1]. Regions are
 *   INDEPENDENT: an ambiguous or missing partner on one side never suppresses the other.</li>
 * <li>Pairing requires same contig, same strand, correct transcript-order side, and
 *   MUTUAL nearest in BOTH directions — with explicit equal-distance TIE DETECTION on both
 *   directions (two candidates at the minimum distance = AMBIGUOUS_PAIRING; first-wins
 *   selection is forbidden; edge-case review requirement 1).</li>
 * <li>Any same-family record INTERSECTING the source record = OVERLAPPING_FLANK for that
 *   family's region: the region emits nothing and pairing NEVER skips past the overlap to
 *   a farther candidate (owner-confirmed disposition). fullITS additionally must not
 *   bridge an overlap-tainted record: a candidate full span containing any record involved
 *   in an OVERLAPPING_FLANK diagnosis is itself rejected.</li>
 * <li>Empty or inverted gaps between NON-intersecting partners = REJECT_EMPTY; no
 *   clipping, no absolute lengths, no invented maximum-distance rule. (Where the computed
 *   gap is inverted BECAUSE the partners intersect, the diagnosis is OVERLAPPING_FLANK —
 *   the cause — not REJECT_*: reconciliation flagged to the fixture author.)</li>
 * <li>Minus-strand intervals keep ascending genomic coordinates and strand '-';
 *   reverse-complementation happens only at sequence extraction (outfasta=).</li>
 * <li>The report= receipt records EVERY pairing decision including rejected alternatives,
 *   so a first-wins regression is visible from the receipt alone (requirement 5).</li>
 * </ul>
 *
 * Input scale note: the GFF loads fully (rRNA rows are few); the optional genome FASTA for
 * outfasta= loads via toReads — fine for fungal genomes; for huge genomes run without
 * outfasta= and extract sequences separately (e.g. cutgff2.sh on the emitted GFF).
 *
 * @author Brian Bushnell
 * @author G11
 * @date 2026-09-09
 */
public class CallITS {

	public static void main(String[] args){
		Timer t=new Timer();
		CallITS x=new CallITS(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public CallITS(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}
		Shared.TRIM_READ_DESCRIPTION=Shared.TRIM_RNAME=true;
		Read.TO_UPPER_CASE=true;
		GffLine.parseAttributes=true;

		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("gff") || a.equals("ingff")){gffIn=b;}
			else if(a.equals("in") || a.equals("fna") || a.equals("ref")){fnaIn=b;}
			else if(a.equals("out") || a.equals("outgff")){gffOut=b;}
			else if(a.equals("report")){reportOut=b;}
			else if(a.equals("outfasta") || a.equals("outfa")){fastaOut=b;}
			else if(a.equals("pattern18s")){pattern18S=b;}
			else if(a.equals("pattern58s") || a.equals("patternr58")){pattern58S=b;}
			else if(a.equals("patternlsu")){patternLSU=b;}
			else if(a.equals("full")){emitFull=Parse.parseBoolean(b);}
			else if(a.equals("individual") || a.equals("its1its2")){emitIndividual=Parse.parseBoolean(b);}
			else if(a.equals("verbose")){verbose=Parse.parseBoolean(b);}
			else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		if(gffIn==null || gffOut==null){
			throw new RuntimeException("CallITS requires gff= and out=. Optional: report=, "
				+"in= + outfasta=, full=t, individual=t, pattern18s/pattern58s/patternlsu.");
		}
		if(fastaOut!=null && fnaIn==null){
			throw new RuntimeException("outfasta= requires in= (the genome fasta).");
		}
	}

	/*--------------------------------------------------------------*/

	public void process(Timer t){
		//Load flank rows. types=null: type filtering is by attribute pattern, since legacy
		//18S rows and generic-RNA rows use different type strings across caller versions.
		ArrayList<GffLine> lines=GffLine.loadGffFile(gffIn, null, false);
		ArrayList<Flank> flanks=new ArrayList<Flank>();
		for(GffLine g : lines){
			if(g.attributes==null){continue;}
			final int fam=(g.attributes.startsWith(pattern18S) ? S18 :
				g.attributes.contains(pattern58S) ? R58 :
				g.attributes.contains(patternLSU) ? LSU : -1);
			if(fam<0){continue;}
			if(g.strand!=GffLine.PLUS && g.strand!=GffLine.MINUS){continue;}//strandless rows unusable
			flanks.add(new Flank(fam, g));
		}
		outstream.println("Loaded "+flanks.size()+" flank rows ("+count(flanks, S18)+" 18S, "
			+count(flanks, R58)+" 5.8S, "+count(flanks, LSU)+" LSU) from "+lines.size()+" GFF lines.");

		//Group by contig+strand; within each group, resolve regions.
		HashMap<String, ArrayList<Flank>> groups=new HashMap<String, ArrayList<Flank>>();
		for(Flank f : flanks){
			String key=f.seqid+"\t"+f.strand;
			ArrayList<Flank> list=groups.get(key);
			if(list==null){groups.put(key, list=new ArrayList<Flank>());}
			list.add(f);
		}

		ArrayList<Region> regions=new ArrayList<Region>();
		ArrayList<String> reportRows=new ArrayList<String>();
		ArrayList<String> keys=new ArrayList<String>(groups.keySet());
		Collections.sort(keys);
		for(String key : keys){
			ArrayList<Flank> group=groups.get(key);
			Collections.sort(group);
			markOverlaps(group, reportRows);
			resolveGroup(group, regions, reportRows);
		}

		writeGff(regions);
		if(reportOut!=null){writeReport(reportRows);}
		if(fastaOut!=null){writeFasta(regions);}

		t.stop();
		long ok=0; for(Region r : regions){if(r.ok()){ok++;}}
		outstream.println("Regions emitted: "+ok+" of "+regions.size()+" evaluated ("
			+count2(regions, "ITS")+" fullITS, "+count2(regions, "ITS1")+" ITS1, "
			+count2(regions, "ITS2")+" ITS2).");
		outstream.println(Tools.timeReadsBasesProcessed(t, flanks.size(), 0, 8));
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/

	/** Diagnoses record intersections, family-aware (derived from the acceptance fixtures):
	 * CROSS-family intersection (e.g. an 18S overlapping the 5.8S pivot) = the
	 * OVERLAPPING_FLANK geometry — both records tainted, affected regions emit nothing,
	 * pairing never skips past. SAME-family intersection (two 5.8S candidates overlapping
	 * each other) = physically impossible as distinct tandem copies — unresolved copy
	 * identity — so BOTH become selfAmbiguous and any pairing involving either (as source
	 * or partner) is AMBIGUOUS_PAIRING; this is what makes the unresolved-copy fixture
	 * reject BOTH candidates even though mutual-nearest alone would pass the closer one. */
	private void markOverlaps(ArrayList<Flank> group, ArrayList<String> report){
		for(int i=0; i<group.size(); i++){
			for(int j=i+1; j<group.size(); j++){
				Flank a=group.get(i), b=group.get(j);
				if(a.start<=b.stop && b.start<=a.stop){//intersecting intervals
					//Both records of every overlap diagnosis carry identity tokens in the
					//receipt (Qiqi rev3 follow-up: the diagnosis row must identify BOTH).
					final String aid=a.brief()+"("+a.id+")", bid=b.brief()+"("+b.id+")";
					if(a.fam==b.fam){
						a.selfAmbiguous=b.selfAmbiguous=true;
						a.ambigWith=(a.ambigWith==null ? bid : a.ambigWith+"|"+bid);
						b.ambigWith=(b.ambigWith==null ? aid : b.ambigWith+"|"+aid);
						report.add(a.seqid+"\t"+a.strand+"\tSAME_FAMILY_OVERLAP(AMBIGUOUS)\t"+aid+"\t"+bid);
					}else{
						a.taintedWith[b.fam]=true;
						b.taintedWith[a.fam]=true;
						report.add(a.seqid+"\t"+a.strand+"\tOVERLAP_DIAGNOSED\t"+aid+"\t"+bid);
					}
				}
			}
		}
	}

	/** Resolves fullITS (per 18S) and ITS1/ITS2 (per 5.8S pivot) for one contig+strand. */
	private void resolveGroup(ArrayList<Flank> group, ArrayList<Region> regions, ArrayList<String> report){
		//fullITS: outer 18S<->LSU pairing, pivot-optional.
		if(emitFull){
			for(Flank f : group){
				if(f.fam!=S18){continue;}
				final boolean plus=(f.strand==GffLine.PLUS);
				//Transcript-downstream side of 18S holds the LSU: genomic-right on plus,
				//genomic-left on minus.
				Pair p=partner(group, f, LSU, plus);
				//Region.left/right are GENOMIC order: on minus the LSU is genomic-left.
				final Flank gl=(p.hit==null || plus ? f : p.hit);
				final Flank gr=(p.hit==null ? null : (plus ? p.hit : f));
				Region r=new Region("ITS", f.seqid, f.strand, gl, gr, null);
				if(p.status!=OK_S){r.status=p.status;}
				else if(f.taintedWith[LSU] || p.hit.taintedWith[S18]){r.status="OVERLAPPING_FLANK";}
				else{
					computeSpan(r, gl, gr, null);
					//No-bridge rule: a full span containing any overlap-involved record is
					//rejected (the affected locus is not a complete path).
					if(r.status.equals(OK_S)){
						for(Flank g2 : group){
							if(g2.anyTaint() && g2.start<=r.stop && r.start<=g2.stop){r.status="OVERLAPPING_FLANK"; break;}
						}
					}
					//Contract decision (Qiqi, 2026-09-09): a DETECTED 5.8S that leaves both
					//spacers empty means the span contains only 5.8S — not an ITS1+5.8S+ITS2
					//product — so it is rejected, not emitted. Distinct from the missing-5.8S
					//rescue, which still emits fullITS from unique outer flanks.
					if(r.status.equals(OK_S)){
						for(Flank g2 : group){
							if(g2.fam==R58 && g2.start>=r.start && g2.stop<=r.stop
								&& g2.start-gl.stop<=1 && gr.start-g2.stop<=1){
								r.status="REJECT_EMPTY"; break;//reason: BOTH_SPACERS_EMPTY
							}
						}
					}
				}
				regions.add(r);
				report.add(r.reportRow(p));
			}
		}
		//ITS1/ITS2: per 5.8S pivot.
		if(emitIndividual){
			for(Flank f : group){
				if(f.fam!=R58){continue;}
				final boolean plus=(f.strand==GffLine.PLUS);
				//ITS1 partner: 18S transcript-upstream (genomic-left on plus, right on minus).
				Pair p1=partner(group, f, S18, !plus);
				final Flank l1=(p1.hit==null ? f : (plus ? p1.hit : f));
				final Flank rr1=(p1.hit==null ? f : (plus ? f : p1.hit));
				Region r1=new Region("ITS1", f.seqid, f.strand, l1, rr1, f);
				if(f.taintedWith[S18]){r1.status="OVERLAPPING_FLANK";}
				else if(p1.status!=OK_S){r1.status=(p1.status=="MISSING" ? "MISSING_18S" : p1.status);}
				else if(p1.hit.taintedWith[R58]){r1.status="OVERLAPPING_FLANK";}
				else{computeSpan(r1, l1, rr1, f);}
				regions.add(r1);
				report.add(r1.reportRow(p1));
				//ITS2 partner: LSU transcript-downstream (genomic-right on plus, left on minus).
				Pair p2=partner(group, f, LSU, plus);
				final Flank l2=(p2.hit==null ? f : (plus ? f : p2.hit));
				final Flank rr2=(p2.hit==null ? f : (plus ? p2.hit : f));
				Region r2=new Region("ITS2", f.seqid, f.strand, l2, rr2, f);
				if(f.taintedWith[LSU]){r2.status="OVERLAPPING_FLANK";}
				else if(p2.status!=OK_S){r2.status=(p2.status=="MISSING" ? "MISSING_LSU" : p2.status);}
				else if(p2.hit.taintedWith[R58]){r2.status="OVERLAPPING_FLANK";}
				else{computeSpan(r2, l2, rr2, f);}
				regions.add(r2);
				report.add(r2.reportRow(p2));
			}
		}
	}

	/**
	 * Finds the partner of family fam for source f on the given genomic side, enforcing:
	 * side correctness (candidate fully right of f when rightSide, else fully left),
	 * equal-distance tie detection (>=2 candidates at min distance = AMBIGUOUS_PAIRING),
	 * and MUTUAL nearest with tie detection in the reverse direction too. Overlapping
	 * candidates are never silently skipped: intersection was pre-diagnosed by
	 * markOverlaps and handled by the caller via taint.
	 */
	private Pair partner(ArrayList<Flank> group, Flank f, int fam, boolean rightSide){
		if(f.selfAmbiguous){
			return new Pair(null, "AMBIGUOUS_PAIRING", "source_unresolved_copy="+f.brief()
				+"("+f.id+")|overlaps="+f.ambigWith);
		}
		Flank best=null; long bd=-1; int ties=0;
		for(Flank c : group){
			if(c.fam!=fam){continue;}
			final long d=(rightSide ? c.start-(long)f.stop : f.start-(long)c.stop);
			if(d<=0){continue;}//not on the required side (intersections handled via taint)
			if(bd<0 || d<bd){bd=d; best=c; ties=1;}
			else if(d==bd){ties++;}
		}
		if(best==null){return new Pair(null, "MISSING");}
		if(ties>1 || best.selfAmbiguous){
			//Enumerate the rejected alternatives (all candidates at the minimum distance,
			//or the unresolved-copy partner) so the receipt shows exactly what tied.
			ByteBuilder alts=new ByteBuilder();
			for(Flank c : group){
				if(c.fam!=fam){continue;}
				final long d=(rightSide ? c.start-(long)f.stop : f.start-(long)c.stop);
				if(d<=0){continue;}
				if(d==bd || (c.selfAmbiguous && c==best)){
					if(alts.length>0){alts.append('|');}
					alts.append(c.brief()).append('(').append(c.id).append(')');
				}
			}
			return new Pair(null, "AMBIGUOUS_PAIRING", "alternatives="+alts);
		}
		//Mutual check with tie detection from the partner's perspective.
		Flank rbest=null; long rbd=-1; int rties=0;
		for(Flank c : group){
			if(c.fam!=f.fam){continue;}
			final long d=(rightSide ? best.start-(long)c.stop : c.start-(long)best.stop);
			if(d<=0){continue;}
			if(rbd<0 || d<rbd){rbd=d; rbest=c; rties=1;}
			else if(d==rbd){rties++;}
		}
		if(rbest!=f || rties>1){
			return new Pair(null, "AMBIGUOUS_PAIRING", "partner="+best.brief()+"("+best.id
				+") mutual_nearest="+(rbest==null ? "none" : rbest.brief()+"("+rbest.id+")")
				+(rties>1 ? " reverse_ties="+rties : ""));
		}
		return new Pair(best, OK_S);
	}

	/** Computes [left.stop+1, right.start-1]; leftF/rightF are the genomic-left and
	 * genomic-right bounding records. Non-intersecting empty gap = REJECT_EMPTY. */
	private void computeSpan(Region r, Flank leftF, Flank rightF, Flank pivot){
		final int lo=leftF.stop+1, hi=rightF.start-1;
		if(hi>=lo){r.start=lo; r.stop=hi; r.status=OK_S;}
		else{r.status="REJECT_EMPTY";}
	}

	/*--------------------------------------------------------------*/

	private void writeGff(ArrayList<Region> regions){
		FileFormat ff=FileFormat.testOutput(gffOut, FileFormat.TEXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("##gff-version 3");
		for(Region r : regions){
			if(!r.ok()){continue;}
			ByteBuilder bb=new ByteBuilder();
			bb.append(r.seqid).tab().append("CallITS").tab().append(r.type).tab();
			bb.append(r.start).tab().append(r.stop).tab().append('.').tab();
			bb.append(r.strand==GffLine.PLUS ? '+' : '-').tab().append('.').tab();
			bb.append("len=").append(r.stop-r.start+1);
			bb.append(";left=").append(r.left.brief()).append('(').append(r.left.id).append(')');
			if(r.right!=null){bb.append(";right=").append(r.right.brief()).append('(').append(r.right.id).append(')');}
			if(r.pivot!=null){bb.append(";pivot=").append(r.pivot.brief()).append('(').append(r.pivot.id).append(')');}
			bsw.println(bb.toBytes());
		}
		errorState|=bsw.poisonAndWait();
	}

	private void writeReport(ArrayList<String> rows){
		FileFormat ff=FileFormat.testOutput(reportOut, FileFormat.TEXT, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		bsw.println("#region\tseqid\tstrand\tstatus\tstart\tstop\tleft(id)\tpartner(id)_or_status\tdetail");
		for(String s : rows){bsw.println(s);}
		errorState|=bsw.poisonAndWait();
	}

	private void writeFasta(ArrayList<Region> regions){
		ArrayList<Read> reads=ReadInputStream.toReads(fnaIn, FileFormat.FA, -1);
		HashMap<String, Read> map=new HashMap<String, Read>();
		for(Read r : reads){map.put(r.id, r);}
		FileFormat ff=FileFormat.testOutput(fastaOut, FileFormat.FA, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		for(Region r : regions){
			if(!r.ok()){continue;}
			Read scaf=map.get(r.seqid);
			if(scaf==null){
				outstream.println("WARNING: contig "+r.seqid+" absent from "+fnaIn+"; sequence skipped.");
				continue;
			}
			if(r.stop>scaf.length()){
				outstream.println("WARNING: "+r.type+" "+r.seqid+":"+r.start+"-"+r.stop
					+" exceeds contig length "+scaf.length()+"; sequence skipped.");
				continue;
			}
			byte[] bases=java.util.Arrays.copyOfRange(scaf.bases, r.start-1, r.stop);
			Read out=new Read(bases, null, r.type+"_"+r.seqid+"_"+r.start+"_"+r.stop+"_"+(r.strand==GffLine.PLUS ? "plus" : "minus"), 1);
			//Reverse-complement is an OUTPUT operation only (contract): coordinates above
			//stay ascending-genomic; the sequence is emitted in transcript orientation.
			if(r.strand==GffLine.MINUS){out.reverseComplement();}
			bsw.println(out);
		}
		errorState|=bsw.poisonAndWait();
	}

	/*--------------------------------------------------------------*/

	private static int count(ArrayList<Flank> list, int fam){
		int n=0; for(Flank f : list){if(f.fam==fam){n++;}} return n;
	}
	private static int count2(ArrayList<Region> list, String type){
		int n=0; for(Region r : list){if(r.ok() && r.type.equals(type)){n++;}} return n;
	}

	/*--------------------------------------------------------------*/

	/** One rRNA flank row. Sort order: start coordinate. */
	static class Flank implements Comparable<Flank> {
		Flank(int fam_, GffLine g){
			fam=fam_; seqid=g.seqid; strand=g.strand; start=g.start; stop=g.stop;
			//Compact provenance token: first attribute field (e.g. "18S,fix=c7_copyA" ->
			//"fix=c7_copyA"; real caller rows -> the leading score/model token). Full
			//attributes stay in the input GFF; the receipt needs identity, not the row.
			String s0=g.attributes;
			int comma=s0.indexOf(',');
			if(comma>=0 && comma+1<s0.length()){s0=s0.substring(comma+1);}
			int comma2=s0.indexOf(',');
			if(comma2>0){s0=s0.substring(0, comma2);}
			id=(s0.length()>48 ? s0.substring(0, 48) : s0);
		}
		String brief(){return FAM_NAMES[fam]+":"+start+"-"+stop;}
		@Override public int compareTo(Flank o){return Integer.compare(start, o.start);}
		final int fam, strand, start, stop;
		final String seqid, id;
		/** Cross-family intersection, per intersecting family: OVERLAPPING_FLANK is
		 * REGION-scoped (an 18S overlapping the pivot kills ITS1 but not ITS2 — the
		 * acceptance fixtures' independent-region requirement). */
		final boolean[] taintedWith=new boolean[3];
		boolean anyTaint(){return taintedWith[0]||taintedWith[1]||taintedWith[2];}
		/** Same-family intersection: unresolved copy identity, pairing AMBIGUOUS. */
		boolean selfAmbiguous=false;
		/** Identity tokens of the same-family record(s) this one intersects. */
		String ambigWith=null;
	}

	/** A pairing attempt result; detail lists competing/tied alternatives for the receipt
	 * (requirement 5: rejected alternatives must be visible, so a first-wins regression
	 * shows in the report). */
	static class Pair {
		Pair(Flank hit_, String status_){this(hit_, status_, null);}
		Pair(Flank hit_, String status_, String detail_){hit=hit_; status=status_; detail=detail_;}
		final Flank hit; final String status; final String detail;
	}

	/** One evaluated region (emitted only when status==OK). */
	static class Region {
		Region(String type_, String seqid_, int strand_, Flank left_, Flank right_, Flank pivot_){
			type=type_; seqid=seqid_; strand=strand_; left=left_; right=right_; pivot=pivot_;
		}
		boolean ok(){return OK_S.equals(status);}
		String reportRow(Pair p){
			return type+"\t"+seqid+"\t"+(strand==GffLine.PLUS ? "+" : "-")+"\t"+status
				+"\t"+(ok() ? start : ".")+"\t"+(ok() ? stop : ".")
				+"\t"+(left!=null ? left.brief()+"("+left.id+")" : ".")
				+"\t"+(p.hit!=null ? p.hit.brief()+"("+p.hit.id+")" : p.status)
				+"\t"+(p.detail!=null ? p.detail : ".");
		}
		final String type, seqid; final int strand;
		final Flank left, right, pivot;
		int start=-1, stop=-1;
		String status="UNSET";
	}

	/*--------------------------------------------------------------*/

	static final int S18=0, R58=1, LSU=2;
	static final String[] FAM_NAMES={"18S", "5.8S", "LSU"};
	static final String OK_S="OK";

	private String gffIn=null, fnaIn=null, gffOut=null, reportOut=null, fastaOut=null;
	private String pattern18S="18S,";
	private String pattern58S="model:r58_";
	private String patternLSU="model:lsu_";
	private boolean emitFull=true;
	private boolean emitIndividual=true;
	private boolean overwrite=true;
	public boolean errorState=false;
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
}
