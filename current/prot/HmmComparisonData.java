package prot;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import structures.ByteBuilder;

/** Prepares reproducible query cohorts for the offline HMM comparison. @author Keqing */
public final class HmmComparisonData {

	public static void main(String[] args){
		try{run(args);}catch(Throwable failure){
			failure.printStackTrace(System.err);
			// A failed producer may own non-daemon ByteStreamWriter threads; never leave a failed CLI hanging.
			System.exit(1);
		}
	}
	private static void run(String[] args){
		if(args.length==1 && args[0].equals("selftest")){selftest(); return;}
		final HashMap<String,String> opt=options(args);
		final String mode=required(opt,"mode");
		if(mode.equals("pool")){pool(opt);}
		else if(mode.equals("sample")){sample(opt);}
		else if(mode.equals("evaluate")){evaluate(opt);}
		else if(mode.equals("regrade")){regrade(opt);}
		else if(mode.equals("matchedprepare")){matchedPrepare(opt);}
		else if(mode.equals("matchedfinish")){matchedFinish(opt);}
		else if(mode.equals("diagnose")){diagnose(opt);}
		else{throw new IllegalArgumentException("Unknown data mode: "+mode);}
	}

	/** Input list has sample, FASTA path; all normalized residues and original IDs remain traceable. */
	private static void pool(HashMap<String,String> opt){
		final ByteStreamWriter fasta=writer(required(opt,"out")), map=writer(required(opt,"map"));
		final HashSet<String> seen=new HashSet<String>();
		int count=0,files=0;
		map.println("query\tsample\toriginal_id\tlength");
		for(String[] row : rows(required(opt,"list"))){
			if(row[0].equals("sample")){continue;}
			if(row.length!=2){throw new IllegalArgumentException("Pool list must have sample and FASTA path");}
			final List<ProteinSequence> proteins=ProteinSearch.readFasta(row[1]);
			for(ProteinSequence p : proteins){
				if(!seen.add(p.id)){throw new IllegalArgumentException("Duplicate source protein: "+p.id);}
				final String id="comparison_g"+count++;
				writeFasta(fasta,id,p.enc);
				map.println(new ByteBuilder().append(id).tab().append(row[0]).tab().append(p.id).tab().append(p.length()));
			}
			files++;
		}
		close(fasta); close(map);
		if(files!=100 || count<20000){throw new IllegalArgumentException("Expected 100 genomes and >=20000 proteins: "+files+", "+count);}
		System.err.println("POOL_PASS files="+files+" proteins="+count);
	}

	/** Samples fixed hash-ranked quotas from current assignments; no-hit is not biological truth. */
	private static void sample(HashMap<String,String> opt){
		final List<ProteinSequence> proteins=ProteinSearch.readFasta(required(opt,"in"));
		final HashMap<String,String> calls=new HashMap<String,String>();
		for(String[] row : rows(required(opt,"calls"))){
			if(row[0].equals("query")){continue;}
			if(calls.put(row[0],row[1])!=null){throw new IllegalArgumentException("Duplicate assignment: "+row[0]);}
		}
		if(calls.size()!=proteins.size()){throw new IllegalArgumentException("Pool and assignment cardinalities differ");}
		final ArrayList<ProteinSequence> yes=new ArrayList<ProteinSequence>(), no=new ArrayList<ProteinSequence>();
		for(ProteinSequence p : proteins){
			final String family=calls.get(p.id);
			if(family==null){throw new IllegalArgumentException("Missing assignment: "+p.id);}
			(family.equals("NA") ? no : yes).add(p);
		}
		final java.util.Comparator<ProteinSequence> order=(a,b)->{
			final int c=Long.compareUnsigned(ReducedAlphabetSeedAssay.stableHash("20261007:"+a.id),
				ReducedAlphabetSeedAssay.stableHash("20261007:"+b.id));
			return c!=0 ? c : a.id.compareTo(b.id);
		};
		yes.sort(order); no.sort(order);
		final int n=Integer.parseInt(opt.getOrDefault("n","10000"));
		if(n<1 || yes.size()<n || no.size()<n){throw new IllegalArgumentException("Insufficient cohorts: "+yes.size()+", "+no.size()+", requested="+n);}
		final ByteStreamWriter out=writer(required(opt,"out")), labels=writer(required(opt,"labels"));
		labels.println("query\tcohort\tproduction_family\tlength");
		for(int i=0; i<n; i++){
			for(int cohort=0; cohort<2; cohort++){
				final ProteinSequence p=(cohort==0 ? yes : no).get(i);
				writeFasta(out,p.id,p.enc);
				labels.println(new ByteBuilder().append(p.id).tab().append(cohort==0 ? "tracked" : "untracked")
					.tab().append(calls.get(p.id)).tab().append(p.length()));
			}
		}
		close(out); close(labels);
		System.err.println("SAMPLE_PASS tracked_pool="+yes.size()+" untracked_pool="+no.size()+" selected="+(2*n));
	}

	/** Grades family assignments against explicit HMM criteria, retaining all curve breakpoints. */
	private static void evaluate(HashMap<String,String> opt){
		final List<ProteinSequence> panel=ProteinSearch.readFasta(required(opt,"in"));
		final HashMap<String,Integer> index=new HashMap<String,Integer>();
		for(int i=0; i<panel.size(); i++){
			if(index.put(panel.get(i).id,i)!=null){throw new IllegalArgumentException("Duplicate panel query");}
		}
		final HashMap<String,String> family=new HashMap<String,String>();
		for(String[] r : rows(required(opt,"manifest"))){
			if(!r[0].equals("active_rank")){family.put("family_"+r[1],r[2]);}
		}
		final String root=required(opt,"root");
		final String[] arms={"hbm","ordinary","blosum"};
		final String[][] prediction=new String[3][panel.size()];
		final double[][] scores=new double[3][panel.size()];
		for(int a=0; a<3; a++){
			java.util.Arrays.fill(scores[a],Double.NEGATIVE_INFINITY);
			for(String[] r : rows(root+"/"+arms[a]+".tsv")){
				if(r[0].equals("query")){continue;}
				final Integer i=index.get(r[0]);
				if(i==null || prediction[a][i]!=null){throw new IllegalArgumentException("Unknown/duplicate result query: "+r[0]);}
				prediction[a][i]=r[1]; scores[a][i]=r[2].equals("NA") ? Double.NEGATIVE_INFINITY : finiteNumber(r[2]);
				if(!r[1].equals("NA") && !Double.isFinite(scores[a][i])){throw new IllegalArgumentException("Assigned result lacks a score: "+r[0]);}
			}
			for(String call : prediction[a]){if(call==null){throw new IllegalArgumentException("Incomplete "+arms[a]+" output");}}
		}
		final ByteStreamWriter summary=writer(required(opt,"out"));
		summary.println("hmm_evalue\tarm\tthreshold\treference_positive\ttp_exact\tfp_exact\tfn_exact\tprecision\trecall\tf1\tany_precision\tany_recall");
		final double[] thresholds={Double.NEGATIVE_INFINITY,Double.NaN,Double.NaN};
		for(double cutoff : new double[]{1e-5,1e-3,1e-10}){
			final String[] truth=new String[panel.size()]; java.util.Arrays.fill(truth,"NA");
			final double[] best=new double[panel.size()]; java.util.Arrays.fill(best,Double.NEGATIVE_INFINITY);
			final ByteFile bf=ByteFile.makeByteFile(required(opt,"domtbl"),false);
			try{
				for(byte[] b=bf.nextLine(); b!=null; b=bf.nextLine()){
					if(b.length==0 || b[0]=='#'){continue;}
					final String[] r=new String(b,StandardCharsets.US_ASCII).trim().split("\\s+");
					if(r.length<22){throw new IllegalArgumentException("Malformed HMMER domain row");}
					final Integer i=index.get(r[0]); final String rep=family.get(r[3]);
					if(i==null || rep==null){throw new IllegalArgumentException("Unknown HMM query/model: "+r[0]+", "+r[3]);}
					final int qlen=Integer.parseInt(r[2]),hlen=Integer.parseInt(r[5]);
					final int hf=Integer.parseInt(r[15]),ht=Integer.parseInt(r[16]),qf=Integer.parseInt(r[17]),qt=Integer.parseInt(r[18]);
					if(qlen!=panel.get(i).length() || hlen<1 || hf<1 || hf>ht || ht>hlen || qf<1 || qf>qt || qt>qlen){
						throw new IllegalArgumentException("Invalid HMM domain coordinates/length: "+r[0]);
					}
					final double covH=(ht-hf+1)/(double)hlen,covQ=(qt-qf+1)/(double)qlen;
					final double sequenceE=finiteNumber(r[6]),domainE=finiteNumber(r[12]),bits=finiteNumber(r[7]);
					if(sequenceE<0 || domainE<0){throw new IllegalArgumentException("Negative HMM E-value");}
					if(sequenceE>cutoff || domainE>cutoff || Math.min(covH,covQ)<0.5){continue;}
					if(bits>best[i] || (bits==best[i] && rep.compareTo(truth[i])<0)){truth[i]=rep; best[i]=bits;}
				}
			}finally{if(bf.close()){throw new RuntimeException("HMM domain read failed");}}
			if(cutoff==1e-5){
				final ByteStreamWriter calls=writer(root+"/hmm.reference.tsv"); calls.println("query\tfamily\tsequence_bits");
				for(int i=0; i<truth.length; i++){calls.println(new ByteBuilder().append(panel.get(i).id).tab().append(truth[i]).tab().append(Double.toString(best[i])));}
				close(calls);
			}
			for(int a=0; a<3; a++){
				if(cutoff==1e-5 && a>0){thresholds[a]=curve(root+"/"+arms[a]+".curve.tsv",truth,prediction[a],scores[a]);}
				final long[] c=counts(truth,prediction[a],scores[a],thresholds[a]);
				final ByteBuilder row=new ByteBuilder().append(Double.toString(cutoff)).tab().append(arms[a]).tab()
					.append(Double.toString(thresholds[a])).tab().append(c[0]).tab().append(c[1]).tab().append(c[2]).tab().append(c[0]-c[1]);
				row.tab().append(ratio(c[1],c[1]+c[2]),8).tab().append(ratio(c[1],c[0]),8).tab()
					.append(ratio(2*c[1],c[0]+c[1]+c[2]),8).tab().append(ratio(c[3],c[1]+c[2]),8).tab().append(ratio(c[3],c[0]),8);
				summary.println(row);
				if(cutoff==1e-5){familyMetrics(root+"/"+arms[a]+".families.tsv",truth,prediction[a],scores[a],thresholds[a]);}
			}
			if(cutoff==1e-5){
				long positive=0,present=0;
				for(String[] r : rows(root+"/ordinary.tsv")){
					if(r[0].equals("query")){continue;}
					final int i=index.get(r[0]);
					if(!truth[i].equals("NA")){
						positive++; if((","+r[9]+",").contains(","+truth[i]+",")){present++;}
					}
				}
				final ByteStreamWriter shortlist=writer(root+"/shortlist.tsv");
				shortlist.println("hmm_positive\ttrue_family_in_unmasked_top50\trecall");
				shortlist.println(new ByteBuilder().append(positive).tab().append(present).tab().append(ratio(present,positive),8));
				close(shortlist);
			}
		}
		close(summary);
		System.err.println("EVALUATION_PASS queries="+panel.size()+" primary=1e-5 coverage=0.5");
	}

	/** Exact score-breakpoint curve; the best F1 is exploratory on this same panel, not a held-out estimate. */
	private static double curve(String path,String[] truth,String[] prediction,double[] score){
		final Integer[] order=new Integer[truth.length];
		long positive=0;
		for(int i=0; i<order.length; i++){order[i]=i; if(!truth[i].equals("NA")){positive++;}}
		java.util.Arrays.sort(order,(a,b)->Double.compare(score[b],score[a]));
		final ByteStreamWriter out=writer(path); out.println("threshold\ttp\tfp\tfn\tprecision\trecall\tf1");
		long tp=0,fp=0; double bestF=-1,bestT=Double.POSITIVE_INFINITY;
		for(int at=0; at<order.length;){
			final double cutoff=score[order[at]];
			if(!Double.isFinite(cutoff)){break;}
			int end=at;
			while(end<order.length && score[order[end]]==cutoff){
				final int i=order[end++];
				if(!prediction[i].equals("NA")){if(prediction[i].equals(truth[i])){tp++;}else{fp++;}}
			}
			final double f=ratio(2*tp,positive+tp+fp);
			out.println(new ByteBuilder().append(cutoff,8).tab().append(tp).tab().append(fp).tab().append(positive-tp)
				.tab().append(ratio(tp,tp+fp),8).tab().append(ratio(tp,positive),8).tab().append(f,8));
			if(f>bestF){bestF=f; bestT=cutoff;}
			at=end;
		}
		close(out); return bestT;
	}

	/** Returns reference positives, exact TP, exact FP and any-hit TP. Wrong-family calls are FP and FN. */
	private static long[] counts(String[] truth,String[] predicted,double[] score,double threshold){
		assert(truth.length==predicted.length && score.length==truth.length) : "Metrics require query-complete matched arrays";
		final long[] out=new long[4];
		for(int i=0; i<truth.length; i++){
			final boolean positive=!truth[i].equals("NA"), accepted=!predicted[i].equals("NA") && score[i]>=threshold;
			if(positive){out[0]++;}
			if(accepted){if(predicted[i].equals(truth[i])){out[1]++;}else{out[2]++;} if(positive){out[3]++;}}
		}
		return out;
	}
	private static double ratio(long n,long d){return d==0 ? 0 : n/(double)d;}

	/** Transposes a frozen query/top50 list into per-model query files; every pair is retained. */
	private static void matchedPrepare(HashMap<String,String> opt){
		final String root=required(opt,"root"),dest=required(opt,"dest");
		final List<ProteinSequence> queries=ProteinSearch.readFasta(root+"/panel.faa");
		final HashMap<String,Integer> qi=new HashMap<String,Integer>(),rank=new HashMap<String,Integer>();
		for(int i=0; i<queries.size(); i++){if(qi.put(queries.get(i).id,i)!=null){throw new IllegalArgumentException("Duplicate query");}}
		final List<String[]> models=rows(required(opt,"manifest"));
		models.remove(0);
		final structures.IntList[] members=new structures.IntList[models.size()];
		for(int i=0; i<models.size(); i++){
			if(Integer.parseInt(models.get(i)[0])!=i || rank.put(models.get(i)[2],i)!=null){throw new IllegalArgumentException("Invalid model rank/representative");}
			members[i]=new structures.IntList();
		}
		final HashSet<String> seen=new HashSet<String>(); long pairs=0;
		for(String[] r : rows(root+"/ordinary.tsv")){
			if(r[0].equals("query")){continue;}
			final Integer q=qi.get(r[0]);
			if(q==null || !seen.add(r[0]) || r.length!=10){throw new IllegalArgumentException("Invalid shortlist query");}
			final String[] refs=r[9].split(",",-1); final HashSet<String> unique=new HashSet<String>();
			if(refs.length!=50){throw new IllegalArgumentException("Not50 candidates: "+r[0]);}
			for(String ref : refs){
				final Integer m=rank.get(ref);
				if(m==null || !unique.add(ref)){throw new IllegalArgumentException("Invalid shortlist family: "+ref);}
				members[m].add(q); pairs++;
			}
		}
		if(seen.size()!=queries.size() || pairs!=50L*queries.size()){throw new IllegalArgumentException("Incomplete pair matrix");}
		final ByteStreamWriter manifest=writer(dest+"/models.tsv"); manifest.println("rank\thmm_name\trep_id\tqueries");
		for(int i=0; i<members.length; i++){
			if(members[i].size==0){continue;}
			final ByteStreamWriter out=writer(dest+"/queries/"+i+".faa");
			for(int j=0; j<members[i].size; j++){final ProteinSequence q=queries.get(members[i].get(j)); writeFasta(out,q.id,q.enc);}
			close(out);
			manifest.println(new ByteBuilder().append(i).tab().append("family_").append(models.get(i)[1]).tab().append(models.get(i)[2]).tab().append(members[i].size));
		}
		close(manifest);
		System.err.println("MATCHED_PAIRS_PASS queries="+queries.size()+" pairs="+pairs);
	}

	/** Ranks the best call per query and admits exactly N, using score then stable query ID for boundary ties. */
	private static void matchedFinish(HashMap<String,String> opt){
		final String root=required(opt,"root"),run=required(opt,"run");
		final int n=Integer.parseInt(opt.getOrDefault("n","10000"));
		final HashMap<String,Integer> index=new HashMap<String,Integer>();
		final ArrayList<String> ids=new ArrayList<String>(),truth=new ArrayList<String>();
		int positives=0;
		for(String[] r : rows(root+"/panel.labels.tsv")){
			if(r[0].equals("query")){continue;}
			if(index.put(r[0],ids.size())!=null){throw new IllegalArgumentException("Duplicate expected query");}
			ids.add(r[0]); truth.add(r[2]); if(!r[2].equals("NA")){positives++;}
		}
		if(n<1 || positives!=n){throw new IllegalArgumentException("Equal-count assay requires N expected positives; found "+positives);}
		final ByteStreamWriter summary=writer(run+"/equal_count.tsv");
		summary.println("method\taccepted\tcorrect_family\twrong_family\taccepted_untracked\tprecision\trecall\tboundary_score\ttied_at_boundary\tties_accepted");
		for(String arm : new String[]{"ordinary","blosum","hbm","hmm"}){
			final Pick[] best=new Pick[ids.size()];
			if(!arm.equals("hmm")){
				for(String[] r : rows(run+"/"+arm+".tsv")){
					if(r[0].equals("query")){continue;}
					final Integer i=index.get(r[0]);
					if(i==null || best[i]!=null || r[1].equals("NA")){throw new IllegalArgumentException("Invalid matched winner row");}
					best[i]=new Pick(i,r[1],finiteNumber(r[2]),0);
				}
				for(Pick p : best){if(p==null){throw new IllegalArgumentException("Missing matched winner");}}
			}else{
				for(String[] model : rows(run+"/models.tsv")){
					if(model[0].equals("rank")){continue;}
					final HashSet<String> requested=new HashSet<String>();
					for(ProteinSequence q : ProteinSearch.readFasta(run+"/queries/"+model[0]+".faa")){requested.add(q.id);}
					if(requested.size()!=Integer.parseInt(model[3])){throw new IllegalArgumentException("Candidate subset count drift");}
					final HashSet<String> reported=new HashSet<String>();
					final ByteFile bf=ByteFile.makeByteFile(run+"/hmm/"+model[0]+".tbl",false);
					try{
						for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
							if(line.length==0 || line[0]=='#'){continue;}
							final String[] r=new String(line,StandardCharsets.US_ASCII).trim().split("\\s+");
							final Integer i=index.get(r[0]);
							if(i==null || !requested.contains(r[0]) || !reported.add(r[0]) || !r[2].equals(model[1])){
								throw new IllegalArgumentException("Unexpected/duplicate HMM pair");
							}
							final double e=finiteNumber(r[4]),bits=finiteNumber(r[5]);
							if(e<0){throw new IllegalArgumentException("Negative HMM E-value");}
							// Printed zero is numeric underflow, not infinite evidence. Saturate at the
							// smallest positive double; finite bits break ties among underflowed E-values.
							final Pick p=new Pick(i,model[2],-Math.log10(Math.max(e,Double.MIN_VALUE)),bits);
							final Pick old=best[i];
							if(old==null || p.score>old.score || (p.score==old.score && (p.secondary>old.secondary ||
								(p.secondary==old.secondary && p.family.compareTo(old.family)<0)))){best[i]=p;}
						}
					}finally{if(bf.close()){throw new RuntimeException("HMM table read failed");}}
				}
			}
			final ArrayList<Pick> ranked=new ArrayList<Pick>();
			for(Pick p : best){
				if(p!=null){
					if(!Double.isFinite(p.score) || !Double.isFinite(p.secondary)){
						throw new IllegalArgumentException("Nonfinite ranked call: "+arm+" "+ids.get(p.query));
					}
					ranked.add(p);
				}
			}
			ranked.sort((a,b)->{
				int c=Double.compare(b.score,a.score); if(c==0){c=Double.compare(b.secondary,a.secondary);}
				return c!=0 ? c : ids.get(a.query).compareTo(ids.get(b.query));
			});
			if(ranked.size()<n){throw new IllegalArgumentException("Fewer than N valid winners: "+arm+" "+ranked.size());}
			final double boundary=ranked.get(n-1).score; int tp=0,wrong=0,negative=0,ties=0,tiesAccepted=0;
			final ByteStreamWriter calls=writer(run+"/"+arm+".accepted.tsv"); calls.println("query\tfamily\tscore\tsecondary_score\texpected_family");
			for(int j=0; j<ranked.size(); j++){
				final Pick p=ranked.get(j); if(p.score==boundary){ties++;}
				if(j>=n){continue;} if(p.score==boundary){tiesAccepted++;}
				final String expected=truth.get(p.query);
				if(expected.equals(p.family)){tp++;}else if(expected.equals("NA")){negative++;}else{wrong++;}
				calls.println(new ByteBuilder().append(ids.get(p.query)).tab().append(p.family).tab().append(Double.toString(p.score)).tab()
					.append(Double.toString(p.secondary)).tab().append(expected));
			}
			close(calls);
			if(tp+wrong+negative!=n){throw new AssertionError("Accepted-call accounting failed");}
			summary.println(new ByteBuilder().append(arm).tab().append(n).tab().append(tp).tab().append(wrong).tab().append(negative).tab()
				.append(tp/(double)n,8).tab().append(tp/(double)positives,8).tab().append(Double.toString(boundary)).tab().append(ties).tab().append(tiesAccepted));
		}
		close(summary); System.err.println("MATCHED_EQUAL_COUNT_PASS accepted_per_method="+n);
	}
	/** Joins logging-only replay events to the original HMM disagreements; no assignments are changed. */
	private static void diagnose(HashMap<String,String> opt){
		final String root=required(opt,"root"),traceRoot=required(opt,"trace"),prefix=required(opt,"prefix");
		final java.util.TreeMap<String,String[]> original=new java.util.TreeMap<String,String[]>();
		final HashMap<String,String> reference=new HashMap<String,String>(),modelMap=new HashMap<String,String>();
		for(String[] r : rows(root+"/hbm.tsv")){if(!r[0].equals("query")){if(original.put(r[0],r)!=null){throw new IllegalArgumentException("Duplicate original query");}}}
		for(String[] r : rows(root+"/hmm.reference.tsv")){if(!r[0].equals("query")){if(reference.put(r[0],r[1])!=null){throw new IllegalArgumentException("Duplicate reference query");}}}
		if(!original.keySet().equals(reference.keySet())){throw new IllegalArgumentException("Original and reference query sets differ");}
		final HashMap<String,HashMap<String,String[]>> traces=new HashMap<String,HashMap<String,String[]>>();
		for(String[] r : rows(traceRoot+"/replay.trace.tsv")){
			if(r.length<2 || !original.containsKey(r[0]) || reference.get(r[0]).equals("NA") || reference.get(r[0]).equals(original.get(r[0])[1])){
				throw new IllegalArgumentException("Trace outside original disagreement set");
			}
			final HashMap<String,String[]> t=traces.computeIfAbsent(r[0],k->new HashMap<String,String[]>());
			if(t.put(r[1],r)!=null){throw new IllegalArgumentException("Duplicate diagnostic event: "+r[0]+" "+r[1]);}
		}
		for(String[] r : rows(required(opt,"manifest"))){if(!r[0].equals("active_rank")){modelMap.put("family_"+r[1],r[2]);}}
		final HashMap<String,double[]> coverage=new HashMap<String,double[]>();
		final ByteFile bf=ByteFile.makeByteFile(root+"/hmm.domtbl",false);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				final String[] r=new String(line,StandardCharsets.US_ASCII).trim().split("\\s+");
				if(r.length<22){throw new IllegalArgumentException("Short HMM domain row");}
				if(!traces.containsKey(r[0]) || !reference.get(r[0]).equals(modelMap.get(r[3]))){continue;}
				final int qlen=Integer.parseInt(r[2]),mlen=Integer.parseInt(r[5]);
				final double qcov=(Integer.parseInt(r[18])-Integer.parseInt(r[17])+1)/(double)qlen;
				final double mcov=(Integer.parseInt(r[16])-Integer.parseInt(r[15])+1)/(double)mlen;
				if(qlen<1 || mlen<1 || qcov<=0 || qcov>1 || mcov<=0 || mcov>1){throw new IllegalArgumentException("Invalid domain coverage");}
				if(finiteNumber(r[6])>1e-5 || finiteNumber(r[12])>1e-5 || qcov<0.5 || mcov<0.5){continue;}
				final double bits=finiteNumber(r[13]); final double[] old=coverage.get(r[0]);
				// One actual qualifying domain, not a union of unrelated domains. Stable coverage tie-break.
				if(old==null || bits>old[3] || (bits==old[3] && (qcov>old[0] || (qcov==old[0] && mcov>old[1])))){
					coverage.put(r[0],new double[]{qcov,mcov,qlen,bits});
				}
			}
		}finally{if(bf.close()){throw new RuntimeException("Domain evidence read failed");}}
		final java.util.TreeMap<String,long[]> counts=new java.util.TreeMap<String,long[]>();
		for(String stage : new String[]{"PRESELECTION","RETRIEVAL_MISS","UNALIGNED_LOOKAHEAD","UNDER_CUTOFF","WRONG_FAMILY"}){counts.put(stage,new long[2]);}
		final java.util.TreeMap<String,structures.DoubleList> distributions=new java.util.TreeMap<String,structures.DoubleList>();
		final ByteStreamWriter out=writer(prefix+".tsv");
		out.println("query\thmm_family\thbm_family\tstage\tfailed_gate\tgate_score\tgate_cutoff\tmargin\ttarget_R\twinner_R\ttarget_hbm_relative\twinner_hbm_relative\tquery_length\thmm_query_coverage\thmm_model_coverage\texploratory_80pct_100aa\tproduction_reason");
		int disagreements=0;
		for(java.util.Map.Entry<String,String[]> entry : original.entrySet()){
			final String q=entry.getKey(),ref=reference.get(q); final String[] call=entry.getValue();
			if(ref.equals("NA") || ref.equals(call[1])){continue;}
			disagreements++;
			final HashMap<String,String[]> t=traces.get(q); final double[] cov=coverage.get(q);
			if(t==null || !t.containsKey("END") || cov==null){throw new IllegalArgumentException("Incomplete diagnostic evidence for "+q);}
			final String[] pre=t.get("PRE"),shortlist=t.get("SHORT"),aln=t.get("ALIGN"),core=t.get("CORE"),hbm=t.get("HBM");
			String stage,gate="NA",value="NA",cutoff="NA";
			if(t.containsKey("NOVALID")){stage="PRESELECTION"; gate="NO_VALID_KMER";}
			else if(pre==null){throw new IllegalArgumentException("Missing target preselection: "+q);}
			else if(!pre[2].equals("true") || !pre[3].equals("true") || !pre[4].equals("true")){
				stage="PRESELECTION"; gate=(!pre[2].equals("true") ? "length;" : "")+(!pre[3].equals("true") ? "kmer_count;" : "")+(!pre[4].equals("true") ? "query_density;" : "");
			}else if(shortlist==null){throw new IllegalArgumentException("Missing surviving shortlist: "+q);}
			else if(Integer.parseInt(shortlist[2])<0){stage="RETRIEVAL_MISS";}
			else if(aln==null){
				if(Integer.parseInt(shortlist[2])<Integer.parseInt(t.get("END")[2])){throw new IllegalArgumentException("Unlogged aligned target: "+q);}
				stage="UNALIGNED_LOOKAHEAD";
			}else if(!t.containsKey("PASS")){
				stage="UNDER_CUTOFF";
				if(!aln[8].equals("true")){gate="acceptance_disabled";}
				else if(finiteNumber(aln[2])<finiteNumber(aln[5])){gate="D55_raw";value=aln[2];cutoff=aln[5];}
				else if(finiteNumber(aln[3])<finiteNumber(aln[6])){gate="identity_percent";value=aln[3];cutoff=aln[6];}
				else if(finiteNumber(aln[4])<finiteNumber(aln[7])){gate="D55_R";value=aln[4];cutoff=aln[7];}
				else if(core!=null && finiteNumber(core[2])<finiteNumber(core[3])){gate="core_overlap";value=core[2];cutoff=core[3];}
				else if(hbm!=null && finiteNumber(hbm[2])<finiteNumber(hbm[3])){gate="HBM_relative";value=hbm[2];cutoff=hbm[3];}
				else{throw new IllegalArgumentException("No measured failing gate explains rejection: "+q);}
			}else{
				if(call[1].equals("NA")){throw new IllegalArgumentException("Accepted target without production winner: "+q);}
				stage="WRONG_FAMILY";
			}
			final boolean strict=cov[0]>=0.8 && cov[1]>=0.8 && cov[2]>=100;
			counts.get(stage)[0]++; if(strict){counts.get(stage)[1]++;}
			String margin="NA";
			if(!value.equals("NA")){
				final double delta=finiteNumber(value)-finiteNumber(cutoff); margin=Double.toString(delta);
				if(!(delta<0)){throw new IllegalArgumentException("Nonnegative failing margin: "+q);}
				addDistribution(distributions,"margin_all_"+gate,delta);
				if(strict){addDistribution(distributions,"margin_strict_"+gate,delta);}
			}
			addDistribution(distributions,"hmm_query_coverage",cov[0]); addDistribution(distributions,"hmm_model_coverage",cov[1]); addDistribution(distributions,"query_length",cov[2]);
			out.println(String.join("\t",q,ref,call[1],stage,gate,value,cutoff,margin,aln==null ? "NA" : aln[4],call[2],
				hbm==null ? "NA" : hbm[2],call[7],Integer.toString((int)cov[2]),Double.toString(cov[0]),Double.toString(cov[1]),Boolean.toString(strict),call[8]));
		}
		close(out);
		if(disagreements!=traces.size()){throw new IllegalArgumentException("Diagnostic row-set mismatch");}
		final ByteStreamWriter summary=writer(prefix+".counts.tsv"); summary.println("stage\tall_rows\tboth_coverage_ge_0.8_length_ge_100");
		for(java.util.Map.Entry<String,long[]> e : counts.entrySet()){summary.println(new ByteBuilder().append(e.getKey()).tab().append(e.getValue()[0]).tab().append(e.getValue()[1]));}
		close(summary);
		final ByteStreamWriter dist=writer(prefix+".deciles.tsv"); dist.println("measurement\tn\tmin\tp10\tp20\tp30\tp40\tp50\tp60\tp70\tp80\tp90\tmax");
		for(java.util.Map.Entry<String,structures.DoubleList> e : distributions.entrySet()){
			final structures.DoubleList values=e.getValue(); java.util.Arrays.sort(values.array,0,values.size);
			final ByteBuilder row=new ByteBuilder().append(e.getKey()).tab().append(values.size);
			for(int p=0; p<=10; p++){row.tab().append(Double.toString(values.array[(int)((values.size-1L)*p/10)]));}
			dist.println(row);
		}
		close(dist); System.err.println("DIAGNOSTIC_TABLE_PASS disagreements="+disagreements);
	}
	private static void addDistribution(java.util.TreeMap<String,structures.DoubleList> all,String name,double value){
		all.computeIfAbsent(name,k->new structures.DoubleList()).add(value);
	}

	private static final class Pick{
		Pick(int query_,String family_,double score_,double secondary_){query=query_;family=family_;score=score_;secondary=secondary_;}
		final int query; final String family; final double score,secondary;
	}

	/** Regrades saved calls against the original production cohorts, without changing old operating points. */
	private static void regrade(HashMap<String,String> opt){
		final String root=required(opt,"root"),prefix=required(opt,"prefix");
		final ArrayList<String> truthList=new ArrayList<String>();
		final HashMap<String,Integer> index=new HashMap<String,Integer>();
		for(String[] r : rows(root+"/panel.labels.tsv")){
			if(r[0].equals("query")){continue;}
			if(r.length!=4 || index.put(r[0],truthList.size())!=null){throw new IllegalArgumentException("Invalid/duplicate panel label");}
			if(!r[1].equals(r[2].equals("NA") ? "untracked" : "tracked")){throw new IllegalArgumentException("Cohort contradicts family: "+r[0]);}
			truthList.add(r[2]);
		}
		if(truthList.isEmpty()){throw new IllegalArgumentException("Empty reference cohort");}
		final String[] truth=truthList.toArray(new String[truthList.size()]);
		final HashMap<String,Double> oldThreshold=new HashMap<String,Double>();
		for(String[] r : rows(root+"/metrics.tsv")){
			if(!r[0].equals("hmm_evalue") && Double.parseDouble(r[0])==1e-5 && !r[1].equals("hbm")){
				oldThreshold.put(r[1],finiteNumber(r[2]));
			}
		}
		final ByteStreamWriter out=writer(prefix+".tsv");
		out.println("arm\toperating_point\tthreshold\tpositive\ttp\tfp\tfn\tprecision\trecall\tf1\tany_precision\tany_recall");
		for(String arm : new String[]{"hbm","ordinary","blosum","hmm"}){
			final String[] prediction=new String[truth.length]; final double[] score=new double[truth.length];
			java.util.Arrays.fill(score,Double.NEGATIVE_INFINITY);
			for(String[] r : rows(root+"/"+(arm.equals("hmm") ? "hmm.reference" : arm)+".tsv")){
				if(r[0].equals("query")){continue;}
				final Integer i=index.get(r[0]);
				if(i==null || prediction[i]!=null){throw new IllegalArgumentException("Unknown/duplicate "+arm+" query: "+r[0]);}
				prediction[i]=r[1];
				if(!r[1].equals("NA")){score[i]=finiteNumber(r[2]);}
				if(arm.equals("hbm") && !r[1].equals(truth[i])){throw new IllegalArgumentException("Current HBM call differs from cohort: "+r[0]);}
			}
			for(String p : prediction){if(p==null){throw new IllegalArgumentException("Incomplete "+arm+" calls");}}
			final boolean tunable=arm.equals("ordinary") || arm.equals("blosum");
			if(tunable && !oldThreshold.containsKey(arm)){throw new IllegalArgumentException("Missing original cutoff for "+arm);}
			final double old=tunable ? oldThreshold.get(arm) : Double.NEGATIVE_INFINITY;
			writeRegrade(out,arm,"original_cutoff",old,counts(truth,prediction,score,old));
			if(tunable){
				final double tuned=curve(prefix+"."+arm+".curve.tsv",truth,prediction,score);
				writeRegrade(out,arm,"HBM_tuned_F1",tuned,counts(truth,prediction,score,tuned));
			}
		}
		close(out);
		System.err.println("HBM_REFERENCE_REGRADE_PASS queries="+truth.length);
	}
	private static void writeRegrade(ByteStreamWriter out,String arm,String policy,double threshold,long[] c){
		out.println(new ByteBuilder().append(arm).tab().append(policy).tab().append(Double.toString(threshold)).tab()
			.append(c[0]).tab().append(c[1]).tab().append(c[2]).tab().append(c[0]-c[1]).tab()
			.append(ratio(c[1],c[1]+c[2]),8).tab().append(ratio(c[1],c[0]),8).tab()
			.append(ratio(2*c[1],c[0]+c[1]+c[2]),8).tab().append(ratio(c[3],c[1]+c[2]),8).tab().append(ratio(c[3],c[0]),8));
	}
	private static double finiteNumber(String text){
		final double value=Double.parseDouble(text);
		if(!Double.isFinite(value)){throw new IllegalArgumentException("Nonfinite measured value: "+text);}
		return value;
	}

	/** Emits explicit per-family confusion and equal-family summaries at the selected operating point. */
	private static void familyMetrics(String path,String[] truth,String[] predicted,double[] score,double threshold){
		final java.util.TreeMap<String,long[]> counts=new java.util.TreeMap<String,long[]>();
		for(int i=0; i<truth.length; i++){
			if(!truth[i].equals("NA")){counts.computeIfAbsent(truth[i],k->new long[3])[0]++;}
			if(!predicted[i].equals("NA") && score[i]>=threshold){
				final long[] c=counts.computeIfAbsent(predicted[i],k->new long[3]); c[1]++;
				if(predicted[i].equals(truth[i])){c[2]++;}
			}
		}
		final ByteStreamWriter out=writer(path); out.println("family\ttruth_positive\ttp\tfp\tfn\tprecision\trecall");
		double macroP=0,macroR=0; int predictedFamilies=0,truthFamilies=0;
		for(java.util.Map.Entry<String,long[]> entry : counts.entrySet()){
			final long[] c=entry.getValue(); final double p=ratio(c[2],c[1]),r=ratio(c[2],c[0]);
			if(c[1]>0){predictedFamilies++; macroP+=p;} if(c[0]>0){truthFamilies++; macroR+=r;}
			out.println(new ByteBuilder().append(entry.getKey()).tab().append(c[0]).tab().append(c[2]).tab().append(c[1]-c[2])
				.tab().append(c[0]-c[2]).tab().append(p,8).tab().append(r,8));
		}
		out.println(new ByteBuilder().append("#macro_precision_over_predicted_families\t").append(predictedFamilies==0 ? 0 : macroP/predictedFamilies,8));
		out.println(new ByteBuilder().append("#macro_recall_over_reference_families\t").append(truthFamilies==0 ? 0 : macroR/truthFamilies,8));
		close(out);
	}

	/** Independently enumerated exact-family confusion: one correct, one wrong family, one false hit, one missed. */
	private static void selftest(){
		final long[] c=counts(new String[]{"a","b","NA","c"},new String[]{"a","a","x","NA"},
			new double[]{1,1,1,Double.NEGATIVE_INFINITY},0);
		if(!java.util.Arrays.equals(c,new long[]{3,1,2,2})){throw new AssertionError("Confusion oracle failed: "+java.util.Arrays.toString(c));}
		System.err.println("HMM_COMPARISON_DATA_SELFTEST_PASS");
	}

	/** Small tabular metadata reader; values are retained for cohort joins and provenance. */
	static List<String[]> rows(String path){
		final ArrayList<String[]> out=new ArrayList<String[]>();
		final ByteFile bf=ByteFile.makeByteFile(path,false);
		try{
			for(byte[] b=bf.nextLine(); b!=null; b=bf.nextLine()){
				if(b.length>0 && b[0]!='#'){out.add(new String(b,StandardCharsets.UTF_8).split("\t",-1));}
			}
		}finally{if(bf.close()){throw new RuntimeException("Read failed: "+path);}}
		return out;
	}

	static HashMap<String,String> options(String[] args){
		final HashMap<String,String> out=new HashMap<String,String>();
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1 || eq==arg.length()-1 || out.put(arg.substring(0,eq),arg.substring(eq+1))!=null){
				throw new IllegalArgumentException("Malformed or duplicate argument: "+arg);
			}
		}
		return out;
	}
	static String required(HashMap<String,String> opt,String key){
		final String value=opt.get(key);
		if(value==null){throw new IllegalArgumentException("Missing "+key+"=");}
		return value;
	}
	static ByteStreamWriter writer(String path){
		if(new java.io.File(path).exists()){throw new IllegalArgumentException("Refusing existing output: "+path);}
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(path,FileFormat.TEXT,null,false,false,false,false));
		out.start(); return out;
	}
	static void close(ByteStreamWriter out){if(out.poisonAndWait()){throw new RuntimeException("Output write failed");}}
	static byte[] ascii(byte[] enc){
		assert(enc!=null && enc.length>0) : "Alignment and FASTA output require a nonempty validated protein";
		final byte[] raw=new byte[enc.length];
		for(int i=0; i<raw.length; i++){raw[i]=enc[i]==Blosum62.X_CODE ? (byte)'X' : AminoAcid.numberToAcid[enc[i]];}
		return raw;
	}
	static void writeFasta(ByteStreamWriter out,String id,byte[] enc){
		out.print(new ByteBuilder().append('>').append(id).nl().append(ascii(enc)).nl());
	}
}
