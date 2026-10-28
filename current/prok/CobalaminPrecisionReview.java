package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.function.Consumer;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins preserved B12 evidence and explicit decisions without editing truth.
 * Unknown full-GA geometry remains pending; optimistic exclusions never add TPs.
 * Feature rows cover caller-overlapping annotations only. Core-touch counts
 * cannot establish absence of annotation beyond the caller interval.
 * @author Brian Bushnell, Ganyu
 */
public final class CobalaminPrecisionReview {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, CobalaminPrecisionReview.class, false);
		String loci=null, context=null, features=null, members=null, models=null, decisions=null, runs=null, neighborhood=null, out=null;
		int expected=-1, events=-1, sensitive=-1, runCount=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("loci")){loci=value;}else if(key.equals("context")){context=value;}
			else if(key.equals("features")){features=value;}else if(key.equals("members")){members=value;}
			else if(key.equals("models")){models=value;}else if(key.equals("decisions")){decisions=value;}
			else if(key.equals("runs")){runs=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("neighborhood")){neighborhood=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else if(key.equals("events")){events=Integer.parseInt(value);}
			else if(key.equals("sensitive")){sensitive=Integer.parseInt(value);}else if(key.equals("runcount")){runCount=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(loci!=null && context!=null && features!=null && members!=null && models!=null && decisions!=null && runs!=null && neighborhood!=null && out!=null
			&& expected>0 && events>=expected && sensitive>=0 && sensitive<=expected && runCount>0,
			"Bind all evidence, explicit decisions and population dimensions");
		Shared.setThreads(1);
		try{
			final CobalaminPrecisionReview x=new CobalaminPrecisionReview();
			x.loci(loci, expected);x.context(context);x.features(features);x.members(members, events);
			x.models(models, sensitive);x.decisions(decisions);x.classify();x.neighborhood(neighborhood);
			x.runs(runs, runCount);x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Evidence inputs        ----------------*/
	/*--------------------------------------------------------------*/

	private void loci(String path, int expected){
		read(path, NcrnaReviewUnion.LOCUS_HEADER, p->{
			final Locus l=new Locus(p);
			require(!byId.contains(l.id), "Duplicate review locus: "+l.id);byId.put(l.id, l);loci.add(l);
		});
		require(loci.size()==expected, "The complete frozen review population is required");
	}
	private void context(String path){
		read(path, NcrnaReviewJoin.CONTEXT_HEADER, p->{
			final Locus l=locus(p);
			require(!l.context && p.termEquals(l.genome, 1) && p.termEquals(l.assembly, 2) && p.termEquals(l.seqid, 3)
				&& p.parseInt(4)==l.start && p.parseInt(5)==l.stop && p.termEquals(l.strand, 6)
				&& p.parseInt(7)==l.stop-l.start+1 && p.termEquals(l.fasta, 24), "Context must preserve exact original prediction identity");
			l.context=true;l.gff=p.parseString(22);l.gffHash=p.parseString(23);l.fastaHash=p.parseString(25);
			require(l.gffHash.matches("[0-9a-f]{20}") && l.fastaHash.matches("[0-9a-f]{20}"), "Source hashes must be sha80");
			l.featuresExpected=p.parseInt(11);l.regionsExpected=p.parseInt(12);
			require(l.featuresExpected>=0 && l.regionsExpected>=0, "Feature counts cannot be negative");
		});
		for(Locus l:loci){require(l.context, "Missing feature context: "+l.id);}
	}
	private void features(String path){
		read(path, NcrnaReviewJoin.OVERLAP_HEADER, p->{
			final Locus l=locus(p);final String kind=p.parseString(8), strand=p.parseString(19);
			require(p.termEquals(l.genome, 1) && p.termEquals(l.assembly, 2) && p.termEquals(l.fasta, 3) && p.termEquals(l.fastaHash, 4)
				&& p.termEquals(l.gff, 5) && p.termEquals(l.gffHash, 6) && p.termEquals(l.seqid, 13), "Feature source binding changed: "+l.id);
			final int start=p.parseInt(16), stop=p.parseInt(17), overlap=p.parseInt(10);
			require(start>0 && start<=l.contigLength && stop>=start && (long)stop-start+1<=l.contigLength
				&& overlap==NcrnaReviewFeatureContext.overlap(l.start, l.stop, start, stop, l.contigLength) && overlap>0
				&& p.parseInt(7)>0 && l.featureLines.add(p.parseString(7)), "Invalid or duplicate original feature geometry");
			require(strand.equals("+") || strand.equals("-") || strand.equals(".") || strand.equals("?"), "Invalid feature strand");
			final String relation=strand.equals(".") || strand.equals("?")?"UNKNOWN":strand.equals(l.strand)?"SAME":"OPPOSITE";
			require(p.termEquals(relation, 11), "Feature strand relation changed");
			final Feature f=new Feature(raw(p), kind, p.parseString(15), p.parseString(21), strand, start, stop);
			l.features.add(f);if(kind.equals("REGION_METADATA")){l.regionsSeen++;}else{l.featuresSeen++;}
			if(!kind.equals("REGION_METADATA") && !kind.equals("PARENT_ANNOTATION") && relation.equals("SAME") && f.b12){l.annotated=true;}
		});
		for(Locus l:loci){require(l.featuresSeen==l.featuresExpected && l.regionsSeen==l.regionsExpected, "Feature rows must reproduce context counts: "+l.id);}
	}
	private void members(String path, int expected){
		read(path, NcrnaReviewJoin.MEMBER_HEADER, p->{
			final String run=p.parseString(0), call=p.parseString(1), key=run+'\t'+call;
			final Locus l=byId.get(p.parseString(2));
			require(l!=null && !run.isEmpty() && !call.isEmpty() && !membership.contains(key) && p.termEquals("APPARENT_FP", 4),
				"Memberships must uniquely preserve original FULL_CM apparent FPs");
			membership.put(key, l);l.members++;
			Population pop=populations.get(run);if(pop==null){pop=new Population();populations.put(run, pop);}pop.fp++;
		});
		require(membership.size()==expected, "Membership denominator changed");
		for(Locus l:loci){require(l.members>0, "Unreferenced review context: "+l.id);}
	}
	private void models(String path, int expectedSensitive){
		final boolean[] seen=new boolean[6];
		read(path, MODEL_HEADER, p->{
			final int model=model(p.parseString(0)), stage=p.termEquals("DEFAULT", 1)?0:p.termEquals("MAX", 1)?1:-1;
			require(model>=0 && stage>=0 && !seen[stage*3+model], "Each of three models needs one DEFAULT and one MAX table");
			seen[stage*3+model]=true;
			final int index=stage*3+model;
			read(p.parseString(2), NcrnaReviewJoin.CM_HEADER, q->{
				final Locus l=locus(q);require(!l.cmSeen[index], "Duplicate per-model locus summary");
				l.cmSeen[index]=true;
				for(int i=0; i<3; i++){l.expected[index][i]=q.parseInt(i+1);}
				require(l.expected[index][0]>=l.expected[index][1] && l.expected[index][1]>=l.expected[index][2] && l.expected[index][2]>=0,
					"CM summary counters must be nested");
				final int full=l.expected[index][2];
				require(q.termEquals(full==0?"NO_FULL_GA_OVERLAPPING_HIT":full==1?"ONE_FULL_GA_OVERLAPPING_HIT":"MULTIPLE_FULL_GA_OVERLAPPING_HITS", 4),
					"CM status disagrees with counts");
			});
			read(p.parseString(3), NcrnaReviewJoin.HIT_HEADER, q->{
				final Locus l=locus(q);require(l.cmSeen[index], "Raw hit missing its model/stage locus summary");
				final Hit h=new Hit(q, l, model, stage);
				require(hitLines.add(index+"\t"+h.line), "Repeated raw CM line for one model/stage");
				l.hits.add(h);l.observed[index][0]++;
				if(h.overlap>0){l.observed[index][1]++;if(stage==1){l.maxOverlap++;}else{l.defaultOverlap++;}}
				if(h.strict){l.observed[index][2]++;}
				if(h.full){l.full++;}if(h.overlap>0 && (h.from!=1 || h.to!=CLEN[model] || h.bits<GA[model])){l.weak++;if(stage==0){l.defaultWeak++;}}
			});
		});
		for(boolean s:seen){require(s, "All six model/stage tables are required");}
		int sensitive=0;
		for(Locus l:loci){
			boolean strictDefault=false;int maxTables=0;
			for(int i=0; i<6; i++){
				if(i<3){require(l.cmSeen[i], "DEFAULT tables must cover every original locus");strictDefault|=l.expected[i][2]>0;}
				else if(l.cmSeen[i]){maxTables++;}
				for(int j=0; j<3; j++){require(l.observed[i][j]==l.expected[i][j], "Raw hits do not reproduce CM summary: "+l.id+" table="+i);}
			}
			require(maxTables==(strictDefault?0:3), "MAX population must be exactly the original unresolved complement");
			l.sensitive=maxTables==3;if(l.sensitive){sensitive++;}
		}
		require(sensitive==expectedSensitive, "Sensitive subset denominator changed");
	}
	private void decisions(String path){
		read(path, DECISION_HEADER, p->{
			final Locus l=locus(p);final String category=p.parseString(1), reason=p.parseString(2);
			require(l.decision==null && (category.equals(REAL) || category.equals(GEOMETRY) || category.equals(CONFLICT)) && !reason.isEmpty(),
				"Explicit decisions must be unique, reasoned supported-locus dispositions");
			l.decision=category;l.reason=reason;
		});
	}
	private void classify(){
		for(Locus l:loci){
			if(l.decision!=null){
				require(l.full>0, "A supported-locus decision requires an overlapping full-span GA hit: "+l.id);
				if(l.decision.equals(GEOMETRY)){require(l.annotated, "Annotated geometry requires an actual same-strand B12 annotation");}
				if(l.decision.equals(REAL)){require(!l.annotated, "An annotated locus cannot silently become a newly unannotated REAL-B12 call");}
				if(l.decision.equals(CONFLICT)){require(coreConflict(l), "Feature-conflict override requires actual same-strand core/feature overlap");}
				l.category=l.decision;
			}else if(l.full>0){l.category=PENDING;l.reason="Overlapping full-GA geometry requires explicit inspection";}
			else if(l.weak>0){l.category=PARTIAL;l.reason="Reported overlapping partial-span or below-GA evidence from either scan";}
			else{
				require(l.sensitive && l.defaultOverlap==0 && l.maxOverlap==0, "Evidence does not fit the declared bins; inspect rather than guess: "+l.id);
				l.category=NONE;l.reason="No overlapping reported DEFAULT or MAX/T0 hit; not proof of biological absence";
			}
		}
	}
	/** Context is descriptive only; it cannot alter the CM disposition or precision. */
	private void neighborhood(String path){
		read(path, CobalaminGenomicContext.HEADER, p->{
			final Locus l=locus(p);require(l.neighborhood==null, "Repeated genomic context would multiply reviewed loci");
			l.overlapClass=p.parseString(1);l.sameHint=bool(p, 8);l.anyHint=bool(p, 13);
			require(!l.overlapClass.isEmpty(), "Every context needs its explicit overlap class");
			for(int col=2; col<=4; col+=2){
				final int bp=p.parseInt(col);final double fraction=p.parseDouble(col+1);
				require(bp>=0 && bp<=l.stop-l.start+1 && Double.isFinite(fraction)
					&& Math.abs(fraction-bp/(double)(l.stop-l.start+1))<1e-7, "CDS overlap fractions must use the original caller length");
			}
			require(bool(p, 6)==(p.parseInt(2)==l.stop-l.start+1), "Containment flag must agree with same-strand CDS overlap");
			for(int col=7; col<=12; col+=5){
				final int ties=p.parseInt(col+4);require(ties>=0 && (ties==0)==p.termEquals("NA", col), "Nearest-CDS absence and tie count disagree");
				if(ties>0){require(p.parseInt(col)>=0 && p.parseInt(col)<=300, "Neighbor must be within the requested downstream range");}
				else{require(!bool(p, col+1), "An absent downstream CDS cannot supply a B12 name hint");}
			}
			l.neighborhood=raw(p);
		});
		for(Locus l:loci){require(l.neighborhood!=null, "Missing neighborhood cannot silently become intergenic: "+l.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Calls and precision    ----------------*/
	/*--------------------------------------------------------------*/

	private void runs(String path, int expected){
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		read(path, RUN_HEADER, p->{
			final Run r=new Run(p.parseString(0));require(!r.name.isEmpty() && names.add(r.name) && populations.contains(r.name),
				"Each requested run must uniquely bind the original review memberships");
			final int count=p.parseInt(4);require(count>0, "Declare the saved FULL_CM union call count");
			read(p.parseString(1), NcrnaReviewUnion.PANEL_CALL_HEADER, q->{
				if(!q.termEquals("FULL_CM", 0) || !q.termEquals("union", 1)){return;}
				final Call c=new Call(q);require(!r.byId.contains(c.id), "Repeated saved call ID");
				c.locus=membership.get(r.name+'\t'+c.id);
				if(c.outcome.equals("APPARENT_FP")){
					final Locus l=c.locus;require(l!=null && q.termEquals(l.genome, 3) && q.termEquals(l.seqid, 4)
						&& q.parseInt(6)==l.start && q.parseInt(7)==l.stop && q.termEquals(l.strand, 8), "FP membership must match exact saved call geometry");
					r.fp++;
				}else{require(c.locus==null, "A reviewed FP cannot become a TP or unscored call without a truth change");}
				r.byId.put(c.id, c);r.calls.add(c);
			});
			require(r.calls.size()==count && r.fp==populations.get(r.name).fp, "Saved call/FP membership denominator changed: "+r.name);
			final String seedCalls=p.parseString(2), seedWindows=p.parseString(3);
			require(seedCalls.equals("NA")==seedWindows.equals("NA"), "Seed call/window tables must be supplied together");
			final String stages=p.parseString(5);
			require(seedCalls.equals("NA")==stages.equals("NA"), "Identity diagnostics and seed windows must be bound together");
			if(!seedCalls.equals("NA")){seedCounts(r, seedCalls, seedWindows);identities(r, stages);}
			runs.add(r);
		});
		require(runs.size()==expected, "Precision run denominator changed");
	}
	private void seedCounts(Run r, String callsPath, String windowsPath){
		read(callsPath, SEED_CALL_HEADER, p->{
			final Call c=r.byId.get(p.parseString(0));
			require(c!=null && !c.seedSeen && p.termEquals(c.family, 1) && p.termEquals(c.outcome, 2), "Seed evidence must join exact saved call/label/outcome");
			c.seedSeen=true;c.windows=p.parseInt(3);c.min=p.parseInt(4);c.max=p.parseInt(5);
			require(c.windows>0 && c.min>0 && c.max>=c.min, "Invalid producing-window support range");
		});
		read(windowsPath, SEED_WINDOW_HEADER, p->{
			final Call c=r.byId.get(p.parseString(0));require(c!=null && p.termEquals(c.family, 3) && p.termEquals(c.outcome, 4)
				&& p.termEquals(c.genome, 1) && p.termEquals(c.seqid, 2) && p.termEquals(c.strand, 5), "Window ownership differs from saved call");
			final int raw=p.parseInt(10), distinct=p.parseInt(11);
			require(raw>=distinct && distinct>0 && c.windowLines.add(p.parseString(6)), "Invalid or duplicate producing window event");
			c.seenWindows++;c.seenMin=Math.min(c.seenMin, distinct);c.seenMax=Math.max(c.seenMax, distinct);
			c.buckets[Math.min(5, distinct)]++;
			r.windows.add(new Window(c, p));
		});
		for(Call c:r.calls){require(c.seedSeen && c.windows==c.seenWindows && c.min==c.seenMin && c.max==c.seenMax,
			"Window events must reproduce every per-call seed range");}
		r.hasSeeds=true;
	}
	/** The current B12 primary path logs bestId for ACCEPT_IDPASS. GFF score is
	 * A+B*length*identity^2, rounded to two decimals; it is not identity.
	 * Match recorded stage lines AND geometry/model, never a threshold or proximity.
	 * Multiple matching values are retained; the selected Orf instance is not logged. */
	private void identities(Run r, String path){
		Collections.sort(r.windows, (a, b)->Long.compare(a.line, b.line));
		for(int i=1; i<r.windows.size(); i++){require(r.windows.get(i-1).line<r.windows.get(i).line,
			"One producing-window event cannot be assigned to multiple emitted calls");}
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long line=1;int next=0;
		try{
			RrnaResourceIO.header(in, NcrnaStageDiagLogSink.HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				line++;if(next==r.windows.size() || line!=r.windows.get(next).line){continue;}
				final Window w=r.windows.get(next++);final Call c=w.call;p.set(row);
				require(p.terms()==21 && p.termEquals("WINDOW", 0) && p.termEquals("ACCEPT_IDPASS", 8),
					"Only source-verified B12 primary IDPASS diagnostics establish bestId identity; inspect any other acceptance path");
				final int strand=c.strand.equals("+")?0:1, length=p.parseInt(4), a=p.parseInt(16), z=p.parseInt(17);
				require(p.termEquals("cobalamin_"+c.family.toLowerCase(Locale.ROOT), 1)
					&& RrnaResourceIO.first(p.parseString(2)).equals(c.seqid) && p.parseInt(3)==strand
					&& p.termEquals(c.model, 12) && p.parseInt(5)==w.pass && p.parseInt(6)==w.start && p.parseInt(7)==w.stop
					&& p.parseInt(9)==w.occurrences && w.stop<length && a>=0 && z>=a && z<length
					&& c.start==(strand==0?a+1:length-z) && c.stop==(strand==0?z+1:length-a),
					"Identity evidence must match the saved call, model and producing window exactly: "+r.name+"/"+c.id);
				w.identity=p.parseFloat(13);require(Float.isFinite(w.identity) && w.identity>=0f && w.identity<=1f,
					"Accepted alignment identity must be a finite fraction, not a score or cutoff");
				w.raw=raw(p);c.identityMin=Math.min(c.identityMin, w.identity);c.identityMax=Math.max(c.identityMax, w.identity);
			}
		}finally{require(!in.close(), "Saved identity diagnostic read failed: "+path);}
		require(next==r.windows.size(), "Every producing window needs its original identity diagnostic");
	}
	private void write(String out){
		RrnaResourceIO.fresh(out, ".loci.tsv", ".hits.tsv", ".features.tsv", ".calls.tsv", ".precision.tsv", ".no_support.tsv", ".no_support_histogram.tsv", ".identity_windows.tsv", ".context_summary.tsv");
		final ByteStreamWriter lr=RrnaResourceIO.writer(out+".loci.tsv"), hr=RrnaResourceIO.writer(out+".hits.tsv"), fr=RrnaResourceIO.writer(out+".features.tsv");
		final ByteStreamWriter cr=RrnaResourceIO.writer(out+".calls.tsv"), pr=RrnaResourceIO.writer(out+".precision.tsv");
		final ByteStreamWriter nr=RrnaResourceIO.writer(out+".no_support.tsv"), hist=RrnaResourceIO.writer(out+".no_support_histogram.tsv");
		final ByteStreamWriter ir=RrnaResourceIO.writer(out+".identity_windows.tsv");
		final ByteStreamWriter cs=RrnaResourceIO.writer(out+".context_summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			lr.println("review_id\tcategory\treason\tdefault_overlap_hits\tmax_overlap_hits\tfull_ga_overlap_hits\tsame_strand_B12_annotation\tdefault_weak_overlap\tdefault_weak_with_no_MAX_overlap");
			hr.println("model_id\tstage\t"+NcrnaReviewJoin.HIT_HEADER+"\tfull_span_ga_overlap\tcall_overlap_rows_touching_core_B12\tcall_overlap_rows_touching_core_other");
			fr.println(NcrnaReviewJoin.OVERLAP_HEADER+"\tB12_dbxref_any_member");
			for(Locus l:loci){
				lr.print(b.clear().append(l.id).tab().append(l.category).tab().append(l.reason).tab().append(l.defaultOverlap).tab().append(l.maxOverlap)
					.tab().append(l.full).tab().append(l.annotated).tab().append(l.defaultWeak>0).tab().append(l.defaultWeak>0 && l.maxOverlap==0).nl());
				for(Feature f:l.features){fr.print(b.clear().append(f.raw).tab().append(f.b12).nl());}
				for(Hit h:l.hits){
					int b12=0, other=0;
					for(Feature f:l.features){if(coreOverlap(l, h, f)){if(f.b12){b12++;}else{other++;}}}
					hr.print(b.clear().append(MODELS[h.model]).tab().append(h.stage==0?"DEFAULT":"MAX").tab().append(h.raw)
						.tab().append(h.full).tab().append(b12).tab().append(other).nl());
				}
			}
			cr.println("run\tcall_id\temitted_family\toriginal_outcome\treview_id\tcategory\toptimistic_excluded\tproducing_windows\tminimum_distinct_keys\tmaximum_distinct_keys\tmodel\tminimum_matching_window_identity\tmaximum_matching_window_identity\tidentity_evidence\tcontext_"+CobalaminGenomicContext.HEADER);
			ir.println("run\tcall_id\tstage_line\tdistinct_seed_keys\t"+NcrnaStageDiagLogSink.HEADER);
			cs.println("run\temitted_family\tcategory\toverlap_class\tsame_strand_nearest_B12_hint_any_tied\tany_strand_nearest_B12_hint_any_tied\tcall_count");
			pr.println("run\toriginal_TP\tall_apparent_FP\tunscored\tREAL_B12\tANNOTATED_GEOMETRY\tFEATURE_CONFLICT\tPARTIAL_BELOW_GA\tNO_SUPPORT\treview_pending\toptimistic_remaining_FP\tpessimistic_precision\toptimistic_precision\tseed_evidence");
			nr.println("run\temitted_family\tNO_SUPPORT_calls\tseed_evidence_calls\tminimum_1_key_calls\tminimum_2plus_key_calls");
			hist.println("run\temitted_family\tdistinct_seed_bucket\tNO_SUPPORT_producing_window_events\tseed_evidence");
			for(Run r:runs){
				for(Window w:r.windows){ir.print(b.clear().append(r.name).tab().append(w.call.id).tab().append(w.line).tab().append(w.distinct).tab().append(w.raw).nl());}
				int tp=0, fp=0, unscored=0, real=0, geometry=0, conflict=0, partial=0, none=0, pending=0;
				final ObjectMap<String,int[]> contexts=new ObjectMap<String,int[]>(String.class, int[].class);final ArrayList<String> contextKeys=new ArrayList<String>();
				final long[][] buckets=new long[3][6];final int[][] counts=new int[3][4];
				for(Call c:r.calls){
					final Locus l=c.locus;final String category=l==null?c.outcome:l.category;
					final boolean excluded=category.equals(REAL) || category.equals(GEOMETRY);
					if(c.outcome.equals("TP")){tp++;}else if(c.outcome.equals("UNSCORED")){unscored++;}else{
						fp++;if(category.equals(REAL)){real++;}else if(category.equals(GEOMETRY)){geometry++;}
						else if(category.equals(CONFLICT)){conflict++;}else if(category.equals(PARTIAL)){partial++;}
						else if(category.equals(NONE)){
							none++;final int f=model(c.family);counts[f][0]++;
							if(c.seedSeen){counts[f][1]++;counts[f][c.min==1?2:3]++;for(int k=1; k<=5; k++){buckets[f][k]+=c.buckets[k];}}
						}else{require(category.equals(PENDING), "Unknown disposition");pending++;}
					}
					if(category.equals(PARTIAL) || category.equals(NONE)){
						final String key=c.family+'\t'+category+'\t'+l.overlapClass+'\t'+l.sameHint+'\t'+l.anyHint;
						int[] count=contexts.get(key);if(count==null){count=new int[1];contexts.put(key, count);contextKeys.add(key);}count[0]++;
					}
					b.clear().append(r.name).tab().append(c.id).tab().append(c.family).tab().append(c.outcome).tab().append(l==null?"NA":l.id)
						.tab().append(category).tab().append(excluded).tab().append(c.seedSeen?Integer.toString(c.windows):"NA")
						.tab().append(c.seedSeen?Integer.toString(c.min):"NA").tab().append(c.seedSeen?Integer.toString(c.max):"NA").tab().append(c.model)
						.tab().append(c.seedSeen?Float.toString(c.identityMin):"NA").tab().append(c.seedSeen?Float.toString(c.identityMax):"NA")
						.tab().append(!c.seedSeen?"UNAVAILABLE":c.identityMin==c.identityMax?"ONE_MATCHING_VALUE":"MULTIPLE_MATCHING_VALUES");
					if(l==null){for(int i=0; i<18; i++){b.append("\tNA");}}else{b.tab().append(l.neighborhood);}cr.print(b.nl());
				}
				Collections.sort(contextKeys);int contextTotal=0;
				for(String key:contextKeys){final int count=contexts.get(key)[0];contextTotal+=count;cs.print(b.clear().append(r.name).tab().append(key).tab().append(count).nl());}
				require(contextTotal==partial+none, "Neighborhood summary must conserve partial and unsupported calls without admitting hints as truth");
				require(tp+fp+unscored==r.calls.size() && fp==real+geometry+conflict+partial+none+pending, "Disposition must conserve original calls and apparent FPs");
				b.clear().append(r.name).tab().append(tp).tab().append(fp).tab().append(unscored).tab().append(real).tab().append(geometry)
					.tab().append(conflict).tab().append(partial).tab().append(none).tab().append(pending).tab();
				if(pending>0){b.append("NA");}else{b.append(fp-real-geometry);}
				b.tab();ratio(b, tp, tp+fp);b.tab();if(pending>0){b.append("NA");}else{ratio(b, tp, tp+fp-real-geometry);}
				pr.print(b.tab().append(r.hasSeeds?"COMPLETE":"PENDING").nl());
				for(int f=0; f<3; f++){
					nr.print(b.clear().append(r.name).tab().append(MODELS[f]).tab().append(counts[f][0]).tab().append(counts[f][1]).tab().append(counts[f][2]).tab().append(counts[f][3]).nl());
					for(int k=1; k<=5; k++){hist.print(b.clear().append(r.name).tab().append(MODELS[f]).tab().append(k==5?"5+":Integer.toString(k))
						.tab().append(r.hasSeeds?Long.toString(buckets[f][k]):"NA").tab().append(r.hasSeeds?"COMPLETE":"PENDING").nl());}
				}
			}
		}finally{RrnaResourceIO.finish(lr, hr, fr, cr, pr, nr, hist, ir, cs);}
		System.err.println("COBALAMIN_PRECISION_REVIEW_WRITTEN loci="+loci.size()+" runs="+runs.size()+"; inspect review_pending; no truth edits or added TPs");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Helpers             ----------------*/
	/*--------------------------------------------------------------*/

	private static boolean coreConflict(Locus l){for(Hit h:l.hits){if(h.full){for(Feature f:l.features){if(!f.b12 && coreOverlap(l, h, f)){return true;}}}}return false;}
	/** Positive overlap in the supplied caller-overlap inventory; false is not
	 * an annotation-absence verdict over core extensions outside the caller. */
	private static boolean coreOverlap(Locus l, Hit h, Feature f){
		return !f.kind.equals("REGION_METADATA") && !f.kind.equals("PARENT_ANNOTATION") && f.strand.equals(h.strand)
			&& NcrnaReviewFeatureContext.overlap(h.start, h.stop, f.start, f.stop, l.contigLength)>0;
	}
	private Locus locus(LineParser1 p){final Locus l=byId.get(p.parseString(0));require(l!=null, "Unknown review ID");return l;}
	private static boolean rfam(String attrs, String model){
		final String db=RrnaResourceIO.optionalField(attrs, "Dbxref");if(db==null){return false;}
		for(int from=0; from<db.length(); ){
			int end=db.indexOf(',', from);if(end<0){end=db.length();}
			if(end-from==5+model.length() && db.regionMatches(true, from, "Rfam:", 0, 5) && db.regionMatches(from+5, model, 0, model.length())){return true;}from=end+1;
		}return false;
	}
	private static int model(String s){for(int i=0; i<MODELS.length; i++){if(MODELS[i].equals(s)){return i;}}return -1;}
	private static void ratio(ByteBuilder b, int numerator, int denominator){if(denominator==0){b.append("NA");}else{b.append(numerator/(double)denominator, 8);}}
	private static byte[] raw(LineParser1 p){return p.line().clone();}
	private static boolean bool(LineParser1 p, int column){require(p.termEquals("true", column) || p.termEquals("false", column), "Malformed boolean");return p.termEquals("true", column);}
	private static void read(String path, String header, Consumer<LineParser1> action){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.UTF_8));final int fields=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==fields, "Changed schema: "+path);action.accept(p);}}
		finally{require(!in.close(), "Input read failed: "+path);}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Locus {
		Locus(LineParser1 p){
			id=p.parseString(0);genome=p.parseString(1);assembly=p.parseString(2);fasta=p.parseString(3);seqid=p.parseString(4);strand=p.parseString(6);
			start=p.parseInt(7);stop=p.parseInt(8);contigLength=p.parseInt(9);
			require(!id.isEmpty() && p.termEquals("cobalamin", 5) && p.termEquals("PREDICTION_UNRESOLVED", 13) && start>0 && stop>=start && stop<=contigLength
				&& (strand.equals("+") || strand.equals("-")), "Invalid original review geometry");
		}
		final String id, genome, assembly, fasta, seqid, strand;final int start, stop, contigLength;
		String gff, gffHash, fastaHash, decision, category, reason, overlapClass;boolean context, annotated, sensitive, sameHint, anyHint;byte[] neighborhood;
		int featuresExpected, regionsExpected, featuresSeen, regionsSeen, members, full, weak, defaultWeak, defaultOverlap, maxOverlap;
		final boolean[] cmSeen=new boolean[6];final int[][] expected=new int[6][3], observed=new int[6][3];
		final ArrayList<Hit> hits=new ArrayList<Hit>();final ArrayList<Feature> features=new ArrayList<Feature>();
		final ObjectSet<String> featureLines=new ObjectSet<String>(String.class);
	}
	private static final class Hit {
		Hit(LineParser1 p, Locus l, int m, int s){
			model=m;stage=s;line=p.parseInt(1);start=p.parseInt(3);stop=p.parseInt(4);strand=p.parseString(5);from=p.parseInt(6);to=p.parseInt(7);bits=p.parseDouble(9);
			// Infernal E-values use exponents; the fast decimal-only parser rejects them.
			p.setBounds(10);final double evalue=Parse.parseDoubleSlow(p.line(), p.a(), p.b());
			require(line>0 && p.termEquals(l.seqid, 2) && start>0 && stop>=start && stop<=l.contigLength && from>=1 && to>=from && to<=CLEN[m]
				&& (strand.equals("+") || strand.equals("-")) && Double.isFinite(bits) && Double.isFinite(evalue) && evalue>=0, "Invalid raw CM geometry/score");
			require(p.termEquals("no", 8) || p.termEquals("5'", 8) || p.termEquals("3'", 8) || p.termEquals("5'&3'", 8), "Unknown CM truncation flag");
			overlap=strand.equals(l.strand)?Math.max(0, Math.min(stop, l.stop)-Math.max(start, l.start)+1):0;
			require(overlap==p.parseInt(12), "Stored sense overlap differs from original prediction");
			full=overlap>0 && from==1 && to==CLEN[m] && bits>=GA[m];
			strict=full && bool(p, 11) && p.termEquals("no", 8) && 2L*overlap>=l.stop-l.start+1 && 2L*overlap>=stop-start+1;
			require(strict==bool(p, 13), "Stored reciprocal-half flag disagrees with raw evidence");
			raw=raw(p);
		}
		final int model, stage, line, start, stop, from, to, overlap;final String strand;final double bits;final boolean full, strict;final byte[] raw;
	}
	private static final class Feature {
		Feature(byte[] r, String k, String t, String attrs, String s, int a, int z){
			raw=r;kind=k;type=t;strand=s;start=a;stop=z;b12=rfam(attrs, MODELS[0]) || rfam(attrs, MODELS[1]) || rfam(attrs, MODELS[2]);
		}
		final byte[] raw;final String kind, type, strand;final int start, stop;final boolean b12;
	}
	private static final class Call {
		Call(LineParser1 p){
			id=p.parseString(2);genome=p.parseString(3);seqid=p.parseString(4);family=p.parseString(5);strand=p.parseString(8);outcome=p.parseString(10);
			model=p.parseString(9);start=p.parseInt(6);stop=p.parseInt(7);
			require(model(family)>=0 && (outcome.equals("TP") || outcome.equals("APPARENT_FP") || outcome.equals("UNSCORED")), "Unknown call subtype/outcome");
			require(start>0 && stop>=start && (strand.equals("+") || strand.equals("-")) && model.startsWith("cobalamin_"+family+'_'),
				"Every saved call must retain exact geometry and its emitting model");
		}
		final String id, genome, seqid, family, strand, outcome, model;final int start, stop;Locus locus;boolean seedSeen;
		float identityMin=Float.POSITIVE_INFINITY, identityMax=Float.NEGATIVE_INFINITY;
		int windows, min, max, seenWindows, seenMin=Integer.MAX_VALUE, seenMax;final long[] buckets=new long[6];
		final ObjectSet<String> windowLines=new ObjectSet<String>(String.class);
	}
	private static final class Window {
		Window(Call c, LineParser1 p){call=c;line=p.parseLong(6);pass=p.parseInt(7);start=p.parseInt(8);stop=p.parseInt(9);occurrences=p.parseInt(10);distinct=p.parseInt(11);
			require(line>1 && (pass==1 || pass==2) && start>=0 && stop>=start, "Producing-window geometry and physical stage line must be valid");}
		final Call call;final long line;final int pass, start, stop, occurrences, distinct;float identity;byte[] raw;
	}
	private static final class Run {
		Run(String n){name=n;}final String name;boolean hasSeeds;int fp;
		final ArrayList<Window> windows=new ArrayList<Window>();
		final ArrayList<Call> calls=new ArrayList<Call>();final ObjectMap<String,Call> byId=new ObjectMap<String,Call>(String.class, Call.class);
	}
	private static final class Population {int fp;}
	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class, Locus.class), membership=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ObjectMap<String,Population> populations=new ObjectMap<String,Population>(String.class, Population.class);
	private final ObjectSet<String> hitLines=new ObjectSet<String>(String.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<Run> runs=new ArrayList<Run>();
	private static final String[] MODELS={"RF00174", "RF01482", "RF01689"};private static final int[] CLEN={189, 149, 81};private static final double[] GA={31.5, 50, 36.5};
	private static final String REAL="REAL-B12", GEOMETRY="ANNOTATED-GEOMETRY", CONFLICT="FEATURE-CONFLICT", PARTIAL="PARTIAL/BELOW-GA", NONE="NO-SUPPORT", PENDING="FULL-GA-REVIEW-REQUIRED";
	public static final String MODEL_HEADER="model_id\tstage\tloci\thits";
	public static final String DECISION_HEADER="review_id\tcategory\treason";
	public static final String RUN_HEADER="run\tcalls\tseed_calls\tseed_windows\texpected_calls\tstages";
	private static final String SEED_CALL_HEADER="call_index\temitted_family\toutcome\tproducing_windows\tminimum_distinct_keys\tmaximum_distinct_keys";
	private static final String SEED_WINDOW_HEADER="call_index\tgenome\tseqid\temitted_family\toutcome\tstrand\tstage_line\tpass\twindow_start0_oriented\twindow_stop0_oriented\tseed_occurrences\tdistinct_seed_keys";
}
