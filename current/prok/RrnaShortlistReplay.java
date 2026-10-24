package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;

/** Reconstructs baseline K9 scores from saved windows, without calling or aligning.
 * True final calls are linked by exact forward coordinates and the GFF model name.
 * Multiple generating windows use the largest required gap, conservatively.
 * Unlinked/ambiguous calls remain in the TP denominator and are never rescued by
 * guessing. FN loci have no observed winning model and are reported separately.
 * @author Raiden
 */
public final class RrnaShortlistReplay {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,RrnaShortlistReplay.class,false);
		String fasta=null,refs=null,stages=null,calls=null,loci=null,prefix=null,assembly=null;
		for(String arg:pp.args){
			final int e=arg.indexOf('=');require(e>0,"Expected flag=value");
			final String k=arg.substring(0,e).toLowerCase(Locale.ROOT),v=arg.substring(e+1);
			if(k.equals("in")){fasta=v;}else if(k.equals("consensus")){refs=v;}
			else if(k.equals("stages")){stages=v;}else if(k.equals("calls")){calls=v;}
			else if(k.equals("loci")){loci=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("assembly")){assembly=v;}else{throw new IllegalArgumentException(k);}
		}
		require(fasta!=null && refs!=null && stages!=null && calls!=null && loci!=null && prefix!=null
			&& assembly!=null,"Require in/consensus/stages/calls/loci/outprefix/assembly");
		final byte[][] library=TrnaConsensusBuilder.loadLibrary(refs);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames.clone();
		new RrnaShortlistReplay(assembly,library,names).run(fasta,stages,calls,loci,prefix);
		Shared.closeStream(pp.outstream);
	}

	public RrnaShortlistReplay(String assembly_,byte[][] library_,String[] names_){
		require(assembly_!=null && !assembly_.isEmpty() && library_!=null && library_.length>1
			&& names_!=null && names_.length==library_.length,"Replay needs an assembly and a multi-model baseline library");
		assembly=assembly_;library=library_;names=names_;
		for(int i=0;i<names.length;i++){
			require(names[i]!=null && names[i].equals(RrnaResourceIO.first(names[i])) && !modelIds.contains(names[i]),"Model names must be unique tokens");
			modelIds.put(names[i],i);
		}
		// Frozen M1 contract: K9, nonadaptive minimum1, whole-library shortlist.
		index=new TrnaKmerIndex(library,9,false,0,0,0,1);
	}

	public void run(String fasta,String stages,String calls,String loci,String prefix){
		RrnaResourceIO.fresh(prefix,".scores.tsv.gz",".links.tsv.gz",".curve.tsv",".summary.tsv");
		loadLoci(loci);loadCalls(calls);loadStages(stages);
		final ByteStreamWriter scores=RrnaResourceIO.writer(prefix+".scores.tsv.gz");
		final ByteBuilder b=new ByteBuilder();scores.println("assembly\tevent\tseqid\tstrand\tpass\twStart0\twStop0\tmodel\tsharedK9\tbestK9\tacceptedModel");
		Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		try{RrnaResourceIO.stream(fasta,r->scoreContig(r,scores,b));}
		finally{RrnaResourceIO.finish(scores);}
		for(Contig c:contigOrder){require(c.seen,"Stage contig absent from pinned FASTA: "+c.name);}
		write(prefix);
		System.err.println("RRNA_SHORTLIST_REPLAY_PASS assembly="+assembly+" windows="+windows
			+" baselineTP="+trueCalls+" linkedTP="+linked+" ambiguous="+ambiguous+" missing="+missing+" alignments=0");
	}

	public void loadLoci(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,LOCI_HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==12 && p.termEquals(assembly,0),"Locus schema/assembly changed");
				if(p.termEquals("TRUTH",2) && p.termEquals("FN",7)){falseNegatives++;}
				if(!p.termEquals("CALL",2) || !p.termEquals("TP",7)){continue;}
				final String seq=p.parseString(1);final char strand=RrnaSeedBenchmarkGrade.strand(p.parseString(6));
				final int start=p.parseInt(4),stop=p.parseInt(5);
				require(start>0 && stop>=start,"Final call coordinates must be forward one-based inclusive");
				final String key=key(seq,strand,start,stop);Call call=trueByKey.get(key);
				if(call==null){call=new Call(seq,strand,start,stop);trueByKey.put(key,call);callOrder.add(call);}
				call.copies++;trueCalls++;
			}
		}finally{require(!in.close(),"Locus read failed");}
	}

	public void loadCalls(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){
			if(row.length==0 || row[0]=='#'){continue;}p.set(row);require(p.terms()==9,"Expected GFF9");
			if(!p.termEquals("rRNA",2) && !p.termEquals("RNA",2)){continue;}
			if(!FAMILY.equals(RrnaSeedBenchmarkGrade.family(row,p,false))){continue;}
			final Call c=trueByKey.get(key(p.parseString(0),RrnaSeedBenchmarkGrade.strand(p.parseString(6)),p.parseInt(3),p.parseInt(4)));
			if(c==null){continue;}c.gffRows++;
			final String model=model(p.parseString(8));require(modelIds.contains(model),"Final call model missing from baseline library: "+model);
			final int id=modelIds.get(model);
			if(c.model>=0 && c.model!=id){c.conflictingModels=true;}else{c.model=id;}
		}}finally{require(!in.close(),"GFF read failed");}
	}

	public void loadStages(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,NcrnaStageDiagLogSink.HEADER,path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==21,"Expected baseline 21-column diagnostic schema");
				if(!p.termEquals("WINDOW",0) || !p.termEquals(FAMILY,1)){continue;}
				if(p.parseInt(10)<0){preIndexWindows++;continue;}
				final String name=RrnaResourceIO.first(p.parseString(2));final int strand=p.parseInt(3),length=p.parseInt(4),pass=p.parseInt(5),start=p.parseInt(6),stop=p.parseInt(7);
				require((strand==0 || strand==1) && length>0 && start>=0 && stop>=start && stop<length
					&& pass==1,"Replay supports the frozen one-pass baseline and exact oriented windows");
				Contig contig=contigs.get(name);
				if(contig==null){contig=new Contig(name,length);contigs.put(name,contig);contigOrder.add(contig);}
				require(contig.length==length,"Conflicting contig lengths in stages");
				final boolean accepted=p.termStartsWith("ACCEPT_",8);final int m=accepted?p.parseInt(11):-1;
				final Window w=new Window(windows++,strand,start,stop,p.parseInt(10),m);contig.windows.add(w);
				if(!accepted){continue;}
				require(m>=0 && m<names.length && p.termEquals(names[m],12),"Accepted model name/index must match the exact baseline library order");
				final int a=p.parseInt(16),z=p.parseInt(17);require(a>=0 && z>=a && z<length,"Accepted call must fit the oriented contig");
				final int forwardStart=(strand==0?a:length-z-1)+1,forwardStop=(strand==0?z:length-a-1)+1;
				final Call c=trueByKey.get(key(name,strand==0?'+':'-',forwardStart,forwardStop));
				if(c!=null && c.model==m){w.call=c;c.windows++;}
			}
		}finally{require(!in.close(),"Stage read failed");}
	}

	void scoreContig(Read read,ByteStreamWriter out,ByteBuilder b){
		final Contig c=contigs.get(RrnaResourceIO.first(read.id));if(c==null){return;}
		require(!c.seen && read.bases!=null && read.length()==c.length,"FASTA contig must bind exactly once with the saved length: "+c.name);c.seen=true;
		for(Window w:c.windows){
			final int a=w.strand==0?w.start:c.length-w.stop-1,z=w.strand==0?w.stop:c.length-w.start-1;
			final byte[] query=Arrays.copyOfRange(read.bases,a,z+1);
			if(w.strand==1){AminoAcid.reverseComplementBasesInPlace(query);}
			final int[] shortlist=index.shortlist(query,names.length);final int best=index.lastMaxShared();
			require(shortlist.length==w.shortlisted,"Reconstructed baseline shortlist size differs; check pinned sequence, model library and config");
			for(int m=0;m<names.length;m++){
				final int score=index.lastSharedCount(m),gap=best-score;
				require(score>=0 && gap>=0 && gap<GAPS,"K9 score outside the production index bound");
				if(score>0){admissionGaps[gap]++;maxGap=Math.max(maxGap,gap);}
				out.print(b.clear().append(assembly).tab().append(w.id).tab().append(c.name).tab().append(w.strand).tab().append(1)
					.tab().append(w.start).tab().append(w.stop).tab().append(names[m]).tab().append(score).tab().append(best).tab().append(m==w.model?1:0).nl());
			}
			if(w.call!=null){
				final int score=index.lastSharedCount(w.model);require(score>0,"An accepted baseline model must have a positive K9 score");
				final int gap=best-score;w.call.minGap=Math.min(w.call.minGap,gap);w.call.maxGap=Math.max(w.call.maxGap,gap);
			}
		}
		// Retain metadata only until its single FASTA contig has been processed.
		c.windows.clear();
	}

	void write(String prefix){
		final ByteStreamWriter links=RrnaResourceIO.writer(prefix+".links.tsv.gz"),curve=RrnaResourceIO.writer(prefix+".curve.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			links.println("assembly\tseqid\tstrand\tstart1\tstop1\tmodel\tTPcopies\tgffRows\twindows\tminGap\tmaxGap\tstatus");
			for(Call c:callOrder){
				final String status;
				if(c.copies!=1 || c.gffRows>1 || c.conflictingModels){status="AMBIGUOUS_FINAL_CALL";ambiguous+=c.copies;}
				else if(c.gffRows!=1 || c.windows==0){status="MISSING_LINK";missing++;}
				else{
					require(c.maxGap>=0 && c.minGap<=c.maxGap,"Every linked call needs actual scored windows");
					status=c.windows==1?"UNIQUE_WINDOW":"MULTIPLE_WINDOWS_CONSERVATIVE";linked++;winnerGaps[c.maxGap]++;
					if(c.windows>1){multipleWindows++;}maxGap=Math.max(maxGap,c.maxGap);
				}
				links.print(b.clear().append(assembly).tab().append(c.seq).tab().append(c.strand).tab().append(c.start).tab().append(c.stop)
					.tab().append(c.model<0?"NA":names[c.model]).tab().append(c.copies).tab().append(c.gffRows).tab().append(c.windows)
					.tab().append(c.maxGap<0?-1:c.minGap).tab().append(c.maxGap).tab().append(status).nl());
			}
			require(linked+ambiguous+missing==trueCalls,"Link classification must conserve every baseline TP");
			curve.println("margin\twindows\tshortlistedModels\tmeanShortlist\tbaselineTP\tlinkedTP\tretainedTP\tretention\tbaselineFN");
			long admitted=0,retained=0;
			for(int gap=0;gap<=maxGap;gap++){
				admitted+=admissionGaps[gap];retained+=winnerGaps[gap];
				if(gap>0 && admissionGaps[gap]==0 && winnerGaps[gap]==0){continue;}
				b.clear().append(gap).tab().append(windows).tab().append(admitted).tab();mean(b,admitted,windows);
				b.tab().append(trueCalls).tab().append(linked).tab().append(retained).tab();mean(b,retained,trueCalls);
				curve.print(b.tab().append(falseNegatives).nl());
			}
			summary.println("assembly\tmappingWindows\tpreIndexWindows\tbaselineTP\tlinkedTP\tmultipleWindowsTP\tambiguousTP\tmissingTP\tbaselineFN\talignments");
			summary.print(b.clear().append(assembly).tab().append(windows).tab().append(preIndexWindows).tab().append(trueCalls).tab().append(linked)
				.tab().append(multipleWindows).tab().append(ambiguous).tab().append(missing).tab().append(falseNegatives).append("\t0\n"));
		}finally{RrnaResourceIO.finish(links,curve,summary);}
	}
	static String model(String attributes){
		String result=null;
		for(int p=0;p<attributes.length();){
			int end=attributes.indexOf(',',p);if(end<0){end=attributes.length();}
			if(attributes.startsWith("model:",p)){require(result==null,"Duplicate GFF model tag");result=attributes.substring(p+6,end);}
			p=end+1;
		}
		require(result!=null && !result.isEmpty(),"Final euk18S call lacks model tag");return result;
	}
	static String key(String seq,char strand,int start,int stop){return seq+'\t'+strand+'\t'+start+'\t'+stop;}
	static void mean(ByteBuilder b,long n,long d){assert(n>=0 && d>=0):"Nonnegative observed counts required";if(d==0){b.append("NA");}else{b.append(n/(double)d,8);}}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Call{
		Call(String s,char st,int a,int z){seq=s;strand=st;start=a;stop=z;}
		final String seq;final char strand;final int start,stop;int copies,gffRows,windows,model=-1,minGap=Integer.MAX_VALUE,maxGap=-1;boolean conflictingModels;
	}
	static final class Window{
		Window(long id_,int s,int a,int z,int n,int m){id=id_;strand=s;start=a;stop=z;shortlisted=n;model=m;}
		final long id;final int strand,start,stop,shortlisted,model;Call call;
	}
	static final class Contig{
		Contig(String n,int len){name=n;length=len;}final String name;final int length;boolean seen;final ArrayList<Window> windows=new ArrayList<Window>();
	}
	public static final String FAMILY="euk18S",LOCI_HEADER="assembly\tseqid\trole\tindex\tstart\tstop\tstrand\toutcome\tpartnerIndex\td5\td3\tcontrolOverlapAnyStrand";
	static final int GAPS=65536;
	final String assembly;final byte[][] library;final String[] names;final TrnaKmerIndex index;
	final ObjectIntMap<String> modelIds=new ObjectIntMap<String>(String.class);
	final ObjectMap<String,Call> trueByKey=new ObjectMap<String,Call>(String.class,Call.class);
	final ObjectMap<String,Contig> contigs=new ObjectMap<String,Contig>(String.class,Contig.class);
	final ArrayList<Call> callOrder=new ArrayList<Call>();final ArrayList<Contig> contigOrder=new ArrayList<Contig>();
	public final long[] admissionGaps=new long[GAPS],winnerGaps=new long[GAPS];
	public long windows,preIndexWindows,trueCalls,falseNegatives,linked,ambiguous,missing,multipleWindows;int maxGap;
}
