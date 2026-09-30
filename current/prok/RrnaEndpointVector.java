package prok;

import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Materializes per-consensus/end rRNA anchor scores and true-end labels.
 * Input partitions are already CM-oriented, assigned, deduplicated and capped.
 * Jitter is sampled from the explicit measured model/end signed-error table.
 * No network training or caller integration.
 * @author Ganyu
 */
public final class RrnaEndpointVector {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaEndpointVector.class,false);
		String in=null,model=null,consensus=null,startTable=null,stopTable=null,histogram=null,prefix=null;
		long expected=-1,seed=-1;int samples=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected explicit key=value");final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(key.equals("in")){in=v;}else if(key.equals("model")){model=v;}else if(key.equals("consensus")){consensus=v;}else if(key.equals("starttable")){startTable=v;}else if(key.equals("stoptable")){stopTable=v;}
			else if(key.equals("histogram")){histogram=v;}else if(key.equals("outprefix")){prefix=v;}else if(key.equals("expected")){expected=Long.parseLong(v);}else if(key.equals("seed")){seed=Long.parseLong(v);}else if(key.equals("samples")){samples=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}}
		require(in!=null&&model!=null&&consensus!=null&&startTable!=null&&stopTable!=null&&histogram!=null&&prefix!=null&&expected>0&&seed>=0&&samples>0,"Declare all inputs, expected loci, samples per locus/end, seed and fresh output prefix");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;ReadWrite.USE_BGZF=true;ReadWrite.ALLOW_NATIVE_BGZF=true;ReadWrite.FORCE_BGZIP=true;
		final RrnaPositionalKmerTable.Table five=RrnaPositionalKmerTable.load(startTable,model),three=RrnaPositionalKmerTable.load(stopTable,model);
		require(five.end.equals("5prime")&&three.end.equals("3prime"),"Bind each endpoint to its own validated k and anchor table");
		new RrnaEndpointVector(model,reference(consensus,model),five,three,distributions(histogram,model),seed,samples).run(in,prefix,expected);
	}

	RrnaEndpointVector(String model_,byte[] consensus_,RrnaPositionalKmerTable.Table five,RrnaPositionalKmerTable.Table three,Distribution[] hist,long seed,int samples_){
		require(consensus_!=null&&consensus_.length>0&&hist.length==2&&samples_>0&&seed>=0,"Explicit reference, two measured histograms and reproducible samples required");
		model=model_;consensus=consensus_;tables=new RrnaPositionalKmerTable.Table[]{five,three};distributions=hist;samples=samples_;
		random=new Random[]{new Random(seed),new Random(seed^0x5DEECE66DL)};rngSeed=seed;
	}

	/** Shared caller/training feature API. Bounds are zero-based inclusive in the
	 * sense-oriented bases. The identity belongs to this SAME raw candidate span.
	 * False means at least one missing k-mer: caller retains raw ends; no NN row.
	 * Scratch/out are worker-local reusable buffers, never shared between threads. */
	public static boolean fill(RrnaPositionalKmerTable.Table table,byte[] bases,int rawStart,int rawStop,int consensusLength,float contigGC,float identity,float[] scratch,float[] out){
		require(table!=null&&bases!=null&&rawStart>=0&&rawStop>=rawStart&&rawStop<bases.length&&consensusLength>0,"Features require a valid inclusive raw candidate span and reference length");
		require(Float.isFinite(contigGC)&&contigGC>=0&&contigGC<=1&&Float.isFinite(identity)&&identity>=0&&identity<=1,"Use measured full-contig GC fraction and finite candidate-span identity, never truth-derived substitutes");
		require(scratch!=null&&scratch.length==SITES&&out!=null&&out.length==INPUTS,"Feature schema is25 scores then3 extras");
		final int rawEnd=table.end.equals("5prime")?rawStart:rawStop;final long anchor=(long)rawEnd+table.offset;
		require(anchor>=Integer.MIN_VALUE&&anchor<=Integer.MAX_VALUE,"Endpoint plus signed anchor offset must not overflow");
		table.scoreWindow(bases,(int)anchor,scratch);boolean usable=true;
		for(int j=0;j<SITES;j++){out[j]=scratch[j];usable&=Float.isFinite(scratch[j]);}
		out[25]=(rawStop-rawStart+1)/(float)consensusLength;out[26]=contigGC;out[27]=identity;return usable;
	}
	/** Class j denotes an END coordinate, even though score j samples its anchor. */
	public static int endpointFromClass(int rawEnd,int classIndex){
		require(classIndex>=0&&classIndex<SITES,"Only25 ordered endpoint offsets can be decoded");return Math.addExact(rawEnd,classIndex-RADIUS);
	}
	static int labelIndex(int truth,int rawEnd){final long j=(long)truth-rawEnd+RADIUS;require(j>=0&&j<SITES,"Clipped jitter must leave the true endpoint reachable");return (int)j;}
	/** Exact candidate slice, shorter sequence as query (consensus wins length ties),
	 * Quantum traceback identity; same helper is available to the caller parity test. */
	public static float candidateIdentity(byte[] bases,int start,int stop,byte[] reference,AlignmentStats stats){
		require(bases!=null&&reference!=null&&reference.length>0&&stats!=null&&start>=0&&stop>=start&&stop<bases.length,"Identity is defined on an actual candidate span");
		final byte[] span=Arrays.copyOfRange(bases,start,stop+1);Tools.toUpperCase(span);byte[] ref=reference;
		for(byte value:reference){if(value>='a'&&value<='z'){ref=reference.clone();Tools.toUpperCase(ref);break;}}
		final byte[] query=span.length<ref.length?span:ref,target=span.length<ref.length?ref:span;
		// QuantumAligner.alignAndTraceStatic seeds bandwidth from stats.rStart;
		// clear prior-candidate state so training/caller order cannot change identity.
		stats.clear();stats.doTrace=true;final float id=QuantumAligner.alignAndTraceStatic(query,target,stats);require(Float.isFinite(id)&&id>=0&&id<=1,"Quantum candidate identity must be finite and fractional");return id;
	}

	void run(String in,String prefix,long expected){
		RrnaResourceIO.fresh(prefix,".5prime.tsv.gz",".3prime.tsv.gz",".audit.tsv.gz",".stats.tsv");
		final ByteStreamWriter[] writers={RrnaResourceIO.writer(prefix+".5prime.tsv.gz"),RrnaResourceIO.writer(prefix+".3prime.tsv.gz")};
		final ByteStreamWriter audit=RrnaResourceIO.writer(prefix+".audit.tsv.gz");final ByteBuilder b=new ByteBuilder();
		final float[] scratch=new float[SITES],features=new float[INPUTS];final ObjectSet<String> ids=new ObjectSet<String>(String.class);final AlignmentStats stats=new AlignmentStats(true);
		try{for(int end=0;end<2;end++){header(writers[end],tables[end],b,rngSeed,samples);}
			audit.println("record_id\tsample\tend\traw_start\traw_stop\ttrue_end\tanchor_center\tlabel_index\tstatus\toutput_row");
			RrnaResourceIO.stream(in,read->{
				RrnaPositionalKmerTable.checkModel(read.id,model);final String id=RrnaResourceIO.first(read.id);require(ids.add(id),"Duplicate locus IDs would duplicate training weight");
				final int left=integerField(read.id,"lflank"),right=integerField(read.id,"rflank");require(read.bases!=null&&left>=0&&right>=0&&(long)left+right<read.length(),"Actual flanks must leave a nonempty sense-oriented gene");
				final int trueStart=left,trueStop=read.length()-right-1;final float gc=floatField(read.id,"gccontig");require(Float.isFinite(gc)&&gc>=0&&gc<=1,"gccontig is source-contig GC fraction, not extracted-window GC");loci++;
				for(int sample=0;sample<samples;sample++){
					final long proposedStart=(long)trueStart+distributions[0].draw(random[0]),proposedStop=(long)trueStop+distributions[1].draw(random[1]);
					final boolean validSpan=proposedStart>=0&&proposedStop>=proposedStart&&proposedStop<read.length();
					final float identity=validSpan?candidateIdentity(read.bases,(int)proposedStart,(int)proposedStop,consensus,stats):Float.NaN;
					for(int end=0;end<2;end++){
						attempts[end]++;final int truth=end==0?trueStart:trueStop;final long rawEnd=end==0?proposedStart:proposedStop;
						final int label=(int)(truth-rawEnd+RADIUS);require(label>=0&&label<SITES,"Histogram clamp keeps labels inside25 classes");
						final boolean usable=validSpan&&fill(tables[end],read.bases,(int)proposedStart,(int)proposedStop,consensus.length,gc,identity,scratch,features);
						long row=-1;String status="KEPT";if(!validSpan){invalidSpan[end]++;status="INVALID_SPAN";}else if(!usable){missing[end]++;status="MISSING_KMER";}
						else{row=kept[end]++;b.clear();for(int j=0;j<INPUTS;j++){if(j>0){b.tab();}b.appendSlow(features[j]);}for(int j=0;j<SITES;j++){b.tab().append(j==label?1:0);}writers[end].print(b.nl());}
						audit.print(b.clear().append(id).tab().append(sample).tab().append(tables[end].end).tab().append(proposedStart).tab().append(proposedStop).tab().append(truth).tab().append(rawEnd+tables[end].offset).tab().append(label).tab().append(status).tab().append(row).nl());
					}
				}
			});require(loci==expected,"Capped input locus count changed; partial outputs are not a completed dataset");
			for(int end=0;end<2;end++){require(attempts[end]==Math.multiplyExact(loci,samples)&&attempts[end]==kept[end]+missing[end]+invalidSpan[end]&&kept[end]>0,"All requested trials must be accounted and each endpoint dataset nonempty");}
		}finally{RrnaResourceIO.finish(writers[0],writers[1],audit);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".stats.tsv");
		try{summary.println("model\tend\tloci\tattempts\tkept\tmissing_kmer\tinvalid_span\thistogram_count\tclipped_low\tclipped_high\tseed\tsamples");
			for(int e=0;e<2;e++){summary.print(b.clear().append(model).tab().append(tables[e].end).tab().append(loci).tab().append(attempts[e]).tab().append(kept[e]).tab().append(missing[e]).tab().append(invalidSpan[e]).tab().append(distributions[e].total).tab().append(distributions[e].low).tab().append(distributions[e].high).tab().append(rngSeed).tab().append(samples).nl());}
		}finally{RrnaResourceIO.finish(summary);}System.err.println("RRNA_ENDPOINT_VECTOR_PASS model="+model+" loci="+loci+" kept="+kept[0]+","+kept[1]);
	}

	static void header(ByteStreamWriter out,RrnaPositionalKmerTable.Table t,ByteBuilder b,long seed,int samples){
		assert(t!=null):"Each dataset is bound to one saved consensus/end table";
		out.println("#dims\t28\t25");out.println("#schema\tRRNA_ENDPOINT_VECTOR_V1");
		out.print(b.clear().append("#model\t").append(t.model).append("\tend\t").append(t.end).append("\tk\t").append(t.k).append("\tanchor_offset\t").append(t.offset).append("\tradius\t12\tseed\t").append(seed).append("\tsamples\t").append(samples).nl());
		out.println("#coordinates\tscore_j=kmer_at(raw_end+anchor_offset+j-12); label_j=end_at(raw_end+j-12); inclusive_sense_coordinates");
		b.clear().append("#columns");for(int j=0;j<SITES;j++){b.tab().append("score_").append(j-RADIUS);}b.append("\tlength_ratio\tcontig_gc\tidentity");for(int j=0;j<SITES;j++){b.tab().append("end_offset_").append(j-RADIUS);}out.print(b.nl());
	}
	static byte[] reference(String path,String model){
		final byte[][] refs=TrnaConsensusBuilder.loadLibrary(path);final String[] names=TrnaConsensusBuilder.lastLoadedNames;byte[] selected=null;
		for(int i=0;i<refs.length;i++){if(RrnaResourceIO.first(names[i]).equals(model)){require(selected==null,"Consensus names must be unique");selected=refs[i];}}
		require(selected!=null&&selected.length>0,"Assigned consensus must exist in the explicit reference library");return selected;
	}
	static Distribution[] distributions(String path,String model){
		final Distribution[] out={new Distribution(),new Distribution()};final ByteFile bf=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(bf,"model\tend\tsigned_error\tcount",path);
			for(byte[] row;(row=bf.nextLine())!=null;){p.set(row);require(p.terms()==4,"Measured histogram schema changed");if(!p.termEquals(model,0)){continue;}
				final int end=p.termEquals("5prime",1)?0:p.termEquals("3prime",1)?1:-1;require(end>=0,"Measured histogram end must be5prime/3prime");out[end].add(p.parseInt(2),p.parseLong(3));}
		}finally{require(!bf.close(),"Measured histogram read failed");}
		require(out[0].total>0&&out[1].total>0,"Missing measured model/end histogram: no pooled, uniform or residual-model fallback");return out;
	}
	static final class Distribution{
		void add(int error,long count){require(count>0&&(!seen||error>last),"Per-model/end signed-error bins must have positive counts and be strictly increasing (no duplicate rows)");seen=true;last=error;total=Math.addExact(total,count);
			final int d=Math.max(-RADIUS,Math.min(RADIUS,error));weights[d+RADIUS]=Math.addExact(weights[d+RADIUS],count);if(error<-RADIUS){low=Math.addExact(low,count);}if(error>RADIUS){high=Math.addExact(high,count);}}
		int draw(Random rng){require(total>0,"Never sample an empty measured histogram");long bits,value;do{bits=rng.nextLong()>>>1;value=bits%total;}while(bits-value+(total-1)<0);long sum=0;
			for(int i=0;i<SITES;i++){sum+=weights[i];if(value<sum){return i-RADIUS;}}throw new IllegalStateException("Histogram draw must land in a conserved count bin");}
		final long[] weights=new long[SITES];long total,low,high;int last;boolean seen;
	}
	static long fieldRange(String header,String key){
		final String token=key+'=';int at=header.indexOf(token);while(at>=0&&at>0&&!Character.isWhitespace(header.charAt(at-1))){at=header.indexOf(token,at+1);}
		require(at>=0,"Missing required extraction field: "+key);final int start=at+token.length();int end=start;while(end<header.length()&&!Character.isWhitespace(header.charAt(end))){end++;}
		require(end>start,"Empty extraction field: "+key);return ((long)start<<32)|(end&0xffffffffL);
	}
	static int integerField(String h,String key){final long r=fieldRange(h,key);return Parse.parseInt(h,(int)(r>>>32),(int)r);}
	static float floatField(String h,String key){final long r=fieldRange(h,key);return Parse.parseFloat(h,(int)(r>>>32),(int)r);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final int RADIUS=12,SITES=25,INPUTS=28;
	final String model;final byte[] consensus;final RrnaPositionalKmerTable.Table[] tables;final Distribution[] distributions;final Random[] random;final long rngSeed;final int samples;
	long loci;final long[] attempts=new long[2],kept=new long[2],missing=new long[2],invalidSpan=new long[2];
}
