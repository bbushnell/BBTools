package prot;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Frozen-pair experiments only; these calls never feed production caches or networks. */
public final class HbmPositionAssay {
	public static void main(String[] args){
		try{
			final HashMap<String,String> o=HmmComparisonData.options(args);
			if("grade".equals(o.get("mode"))){grade(o);}else if("gradeaccepted".equals(o.get("mode"))){gradeAccepted(o);}else{run(o);}
		}catch(Throwable failure){failure.printStackTrace();System.exit(1);}
	}
	private static void run(HashMap<String,String> o) throws Exception{
		final long start=System.nanoTime();
		final String resources=HmmComparisonData.required(o,"resources"),prefix=HmmComparisonData.required(o,"out");
		final String kind=o.getOrDefault("kind","logodds"),path=o.getOrDefault("path","both");
		if(!path.equals("both") && !path.equals("profile") && !path.equals("blosum")){throw new IllegalArgumentException("Invalid path mode");}
		final boolean validation=parse.Parse.parseBoolean(o.getOrDefault("validate","true"));
		final int threads=Integer.parseInt(o.getOrDefault("t","1"));
		if(threads<1){throw new IllegalArgumentException("Positive worker count required");}
		final List<ProteinSequence> queries=ProteinSearch.readFasta(HmmComparisonData.required(o,"in"));
		final List<ProteinSequence> refs=ProteinSearch.readFasta(o.getOrDefault("refs",resources+"/consensus_reps_round2.faa.gz"));
		final HashMap<String,Integer> ranks=new HashMap<String,Integer>();
		final ArrayList<String> roster=new ArrayList<String>(); final HashMap<String,byte[]> sequence=new HashMap<String,byte[]>();
		for(int i=0; i<refs.size(); i++){
			final ProteinSequence p=refs.get(i);
			if(ranks.put(p.id,i)!=null){throw new IllegalArgumentException("Duplicate model");}
			roster.add(p.id);sequence.put(p.id,p.enc);
		}
		final HashMap<String,int[]> pairs=new HashMap<String,int[]>();
		for(String[] row : HmmComparisonData.rows(HmmComparisonData.required(o,"pairs"))){
			if(row[0].equals("query")){continue;}
			if(row.length!=10){throw new IllegalArgumentException("Frozen shortlist needs ten fields");}
			final String[] names=row[9].split(",",-1);final int[] list=new int[names.length];
			if(names.length!=50){throw new IllegalArgumentException("Expected full top50 for "+row[0]);}
			final HashSet<String> seen=new HashSet<String>();
			for(int j=0; j<names.length; j++){
				final Integer rank=ranks.get(names[j]);
				if(rank==null || !seen.add(names[j])){throw new IllegalArgumentException("Unknown/duplicate model in shortlist");}
				list[j]=rank;
			}
			if(pairs.put(row[0],list)!=null){throw new IllegalArgumentException("Duplicate query shortlist");}
		}
		final HashSet<String> uniqueQueries=new HashSet<String>();
		for(ProteinSequence q : queries){if(!uniqueQueries.add(q.id) || !pairs.containsKey(q.id)){throw new IllegalArgumentException("Duplicate/unlisted query: "+q.id);}}
		if(queries.size()!=pairs.size()){throw new IllegalArgumentException("Frozen query sets differ");}
		final HbmPositionModel[] models;
		if(path.equals("blosum")){models=null;}
		else{
			final HbmBundleLoader.Loaded hbm=HbmBundleLoader.load(java.nio.file.Paths.get(o.getOrDefault("hbm",resources+"/magqc_hbm_v1.rare01.hbmt.gz")),roster,
				id->sequence.get(id),HbmBundleLoader.loadSemanticProvenance(o.getOrDefault("provenance",resources+"/PROVENANCE_MANIFEST.tsv.gz")));
			final double beta=Double.parseDouble(o.getOrDefault("beta","0.1")),floor=Double.parseDouble(o.getOrDefault("clipmin","-4"));
			final boolean clip=parse.Parse.parseBoolean(o.getOrDefault("clip","false"));
			if(o.containsKey("background")){
				final List<String[]> rows=HmmComparisonData.rows(o.get("background"));
				if(rows.size()!=21 || !String.join("\t",rows.get(0)).equals("residue_code\tprobability")){throw new IllegalArgumentException("Frozen background needs its header and exactly twenty residue rows");}
				final double[] background=new double[20];
				for(int a=0; a<20; a++){
					final String[] row=rows.get(a+1);
					if(row.length!=2 || Integer.parseInt(row[0])!=a){throw new IllegalArgumentException("Frozen background residue order differs");}
					background[a]=Double.parseDouble(row[1]);
				}
				models=hbm.positionModels(kind,beta,clip,floor,background);
			}else{models=hbm.positionModels(kind,beta,clip,floor);}
		}
		final long loaded=System.nanoTime();
		final Winner[][] winners=new Winner[queries.size()][2];
		final long[][] work=new long[threads][5];
		final AtomicInteger next=new AtomicInteger();final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final Thread[] workers=new Thread[threads];
		for(int w=0; w<threads; w++){
			final int worker=w;
			workers[w]=new Thread(()->{
				try{
					for(int i=next.getAndIncrement(); i<queries.size() && failure.get()==null; i=next.getAndIncrement()){
						final ProteinSequence q=queries.get(i);
						for(int rank : pairs.get(q.id)){
							final ProteinSequence ref=refs.get(rank);int fixed=Integer.MIN_VALUE;
							if(!path.equals("profile")){
								final long before=System.nanoTime();
								final AAAlignment old=GlocalAminoLinear.align(q.enc,ref.enc,!path.equals("blosum"));
								fixed=path.equals("blosum") ? old.rawScore*64 : models[rank].scorePath(q.enc,old.tStart,old.match);
								work[worker][0]+=System.nanoTime()-before;work[worker][2]++;
								winners[i][0]=better(winners[i][0],new Winner(ref.id,fixed));
							}
							if(models!=null){
								final long before=System.nanoTime();
								final HbmPositionModel.Result result=models[rank].align(q.enc,validation);
								work[worker][1]+=System.nanoTime()-before;work[worker][3]++;
								if(validation && models[rank].scorePath(q.enc,result.start,result.path)!=result.score){throw new AssertionError("Profile traceback score differs from DP");}
								if(!path.equals("profile") && result.score<fixed){throw new AssertionError("Profile optimum is worse than feasible BLOSUM path");}
								if(!path.equals("profile") && result.score>fixed){work[worker][4]++;}
								winners[i][1]=better(winners[i][1],new Winner(ref.id,result.score));
							}
						}
						if(i%5000==0){System.err.println("POSITION_PROGRESS "+i+"/"+queries.size());}
					}
				}catch(Throwable ex){failure.compareAndSet(null,ex);}
			}); workers[w].start();
		}
		for(Thread worker : workers){worker.join();}
		if(failure.get()!=null){throw new RuntimeException("Position scoring worker failed",failure.get());}
		final long scored=System.nanoTime();
		final ByteStreamWriter out=HmmComparisonData.writer(prefix+".tsv");
		out.println("query\tfixed_family\tfixed_score64\tprofile_family\tprofile_score64\tquery_length");
		for(int i=0; i<queries.size(); i++){
			final ByteBuilder row=new ByteBuilder().append(queries.get(i).id);
			for(Winner winner : winners[i]){
				row.tab().append(winner==null ? "NA" : winner.family).tab().append(winner==null ? "NA" : Integer.toString(winner.score));
			}
			out.println(row.tab().append(queries.get(i).length()));
		}
		HmmComparisonData.close(out);
		final long[] counts=new long[5];for(long[] row : work){for(int i=0; i<5; i++){counts[i]+=row[i];}}
		final long expected=50L*queries.size();
		if((!path.equals("profile") && counts[2]!=expected) || (models!=null && counts[3]!=expected)){throw new AssertionError("Position assay did not score every frozen pair");}
		final ByteStreamWriter timing=HmmComparisonData.writer(prefix+".time.tsv");
		timing.println("queries\tthreads\tload_s\tcompute_s\tfixed_worker_elapsed_s\tprofile_worker_elapsed_s\tfixed_pairs\tprofile_pairs\tstrict_profile_improvements");
		timing.println(new ByteBuilder().append(queries.size()).tab().append(threads).tab().append((loaded-start)/1e9,6).tab().append((scored-loaded)/1e9,6)
			.tab().append(counts[0]/1e9,6).tab().append(counts[1]/1e9,6).tab().append(counts[2]).tab().append(counts[3]).tab().append(counts[4]));
		HmmComparisonData.close(timing);System.err.println("POSITION_RUN_PASS pairs="+expected);
	}
	private static Winner better(Winner old,Winner next){
		return old==null || next.score>old.score || (next.score==old.score && next.family.compareTo(old.family)<0) ? next : old;
	}
	private static final class Winner{
		Winner(String f,int s){family=f;score=s;}
		final String family;final int score;
	}
	/** Exact accepted-count comparison against both operational reference assignments. */
	private static void grade(HashMap<String,String> o){
		final String root=HmmComparisonData.required(o,"root"),input=HmmComparisonData.required(o,"in");
		final int accepted=Integer.parseInt(o.getOrDefault("n","10000"));
		final HashMap<String,String> hbm=new HashMap<String,String>(),hmm=new HashMap<String,String>();
		for(String[] r : HmmComparisonData.rows(root+"/panel.labels.tsv")){if(!r[0].equals("query")){if(hbm.put(r[0],r[2])!=null){throw new IllegalArgumentException("Duplicate HBM reference query");}}}
		for(String[] r : HmmComparisonData.rows(root+"/hmm.reference.tsv")){if(!r[0].equals("query")){if(hmm.put(r[0],r[1])!=null){throw new IllegalArgumentException("Duplicate HMM reference query");}}}
		if(!hbm.keySet().equals(hmm.keySet())){throw new IllegalArgumentException("Reference query sets differ");}
		final List<String[]> rows=HmmComparisonData.rows(input);rows.remove(0);
		final HashSet<String> seen=new HashSet<String>();
		for(String[] row : rows){if(row.length!=6 || !hbm.containsKey(row[0]) || !seen.add(row[0]) || Integer.parseInt(row[5])<1){throw new IllegalArgumentException("Invalid position winner row");}}
		if(seen.size()!=hbm.size() || accepted<1 || accepted>seen.size()){throw new IllegalArgumentException("Incomplete panel or invalid accepted count");}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(o,"out"));
		out.println("path\tranking\treference\taccepted\tpositives\ttp\twrong_family\taccepted_reference_negative\tprecision\trecall\tf1\tboundary\tboundary_ties\tties_accepted");
		for(int which=0; which<2; which++){
			final int familyColumn=1+2*which,scoreColumn=familyColumn+1;
			if(rows.get(0)[familyColumn].equals("NA")){continue;}
			for(boolean mean : new boolean[]{false,true}){
				final ArrayList<Ranked> order=new ArrayList<Ranked>();
				for(String[] row : rows){
					if(row[familyColumn].equals("NA")){throw new IllegalArgumentException("Incomplete scoring path");}
					final double score=Integer.parseInt(row[scoreColumn])/(64.0*(mean ? Integer.parseInt(row[5]) : 1));
					order.add(new Ranked(row[0],row[familyColumn],score));
				}
				order.sort((a,b)->{final int c=Double.compare(b.score,a.score);return c!=0 ? c : a.query.compareTo(b.query);});
				final double boundary=order.get(accepted-1).score;int ties=0,tiesAccepted=0;
				for(int i=0; i<order.size(); i++){if(order.get(i).score==boundary){ties++;if(i<accepted){tiesAccepted++;}}}
				for(String reference : new String[]{"production_hbm","hmm"}){
					final HashMap<String,String> truth=reference.equals("hmm") ? hmm : hbm;
					int positive=0,tp=0,wrong=0,negative=0;
					for(String family : truth.values()){if(!family.equals("NA")){positive++;}}
					for(int i=0; i<accepted; i++){
						final Ranked call=order.get(i);final String expected=truth.get(call.query);
						if(expected.equals(call.family)){tp++;}else if(expected.equals("NA")){negative++;}else{wrong++;}
					}
					if(tp+wrong+negative!=accepted){throw new AssertionError("Accepted-call conservation failed");}
					out.println(new ByteBuilder().append(which==0 ? "fixed" : "profile").tab().append(mean ? "per_residue" : "total").tab().append(reference)
						.tab().append(accepted).tab().append(positive).tab().append(tp).tab().append(wrong).tab().append(negative)
						.tab().append(tp/(double)accepted,8).tab().append(positive==0 ? 0 : tp/(double)positive,8)
						.tab().append(2.0*tp/(accepted+positive),8).tab().append(boundary,8).tab().append(ties).tab().append(tiesAccepted));
				}
			}
		}
		HmmComparisonData.close(out);System.err.println("POSITION_GRADE_PASS queries="+rows.size());
	}
	private static final class Ranked{
		Ranked(String q,String f,double s){query=q;family=f;score=s;}
		final String query,family;final double score;
	}
	/** Regrades previously frozen accepted calls, preserving their original ranking and cutoff. */
	private static void gradeAccepted(HashMap<String,String> o){
		final String root=HmmComparisonData.required(o,"root");
		final int n=Integer.parseInt(o.getOrDefault("n","10000"));
		final List<String[]> calls=HmmComparisonData.rows(HmmComparisonData.required(o,"in"));calls.remove(0);
		if(calls.size()!=n){throw new IllegalArgumentException("Frozen accepted-call count differs from N");}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(o,"out"));
		out.println("reference\taccepted\tpositives\ttp\twrong_family\taccepted_reference_negative\tprecision\trecall\tf1");
		for(String reference : new String[]{"production_hbm","hmm"}){
			final HashMap<String,String> truth=new HashMap<String,String>();int positive=0;
			for(String[] row : HmmComparisonData.rows(root+(reference.equals("hmm") ? "/hmm.reference.tsv" : "/panel.labels.tsv"))){
				if(row[0].equals("query")){continue;}
				final String family=row[reference.equals("hmm") ? 1 : 2];
				if(truth.put(row[0],family)!=null){throw new IllegalArgumentException("Duplicate reference query");}
				if(!family.equals("NA")){positive++;}
			}
			int tp=0,wrong=0,negative=0;final HashSet<String> seen=new HashSet<String>();
			for(String[] call : calls){
				final String expected=truth.get(call[0]);
				if(expected==null || !seen.add(call[0]) || call[1].equals("NA")){throw new IllegalArgumentException("Invalid frozen accepted call");}
				if(expected.equals(call[1])){tp++;}else if(expected.equals("NA")){negative++;}else{wrong++;}
			}
			if(tp+wrong+negative!=n){throw new AssertionError("Frozen-call accounting failed");}
			out.println(new ByteBuilder().append(reference).tab().append(n).tab().append(positive).tab().append(tp).tab().append(wrong).tab().append(negative)
				.tab().append(tp/(double)n,8).tab().append(positive==0 ? 0 : tp/(double)positive,8).tab().append(2.0*tp/(n+positive),8));
		}
		HmmComparisonData.close(out);
	}
}
