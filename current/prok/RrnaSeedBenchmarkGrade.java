package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.LongHashSet;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader.MatchResult;
import prok.NcrnaGffFamilyGrader.MatchedPair;
import structures.ByteBuilder;

/** Per-assembly seed-arm grading of final integrated rRNA calls against CM truth.
 * Keeps euk5S calls on prok controls as FP, never prok5S true positives.
 * Reuses the canonical maximum-cardinality reciprocal-half matcher per contig.
 * @author Raiden
 */
public final class RrnaSeedBenchmarkGrade {

	public static void main(String[] args){
		String calls=null,truth=null,prefix=null,arm=null,clade=null,assembly=null,family="euk5S";boolean negative=false;
		for(String arg:args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value: "+arg);final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("calls")){calls=v;}else if(k.equals("truth")){truth=v;}else if(k.equals("outprefix")){prefix=v;}
			else if(k.equals("arm")){arm=v;}else if(k.equals("clade")){clade=v;}else if(k.equals("assembly")){assembly=v;}
			else if(k.equals("family")){family=v;}else if(k.equals("negative")){negative=Parse.parseBoolean(v);}else{throw new IllegalArgumentException(k);}}
		require(calls!=null&&truth!=null&&prefix!=null&&arm!=null&&clade!=null&&assembly!=null,"Required calls/truth/outprefix/arm/clade/assembly");
		require(family.equals("euk5S")||family.equals("prok5S"),"Only explicit euk5S/prok5S grading is supported");
		require(!negative||family.equals("euk5S"),"Negative mode is the euk5S-on-prok cross-family control");
		final RrnaSeedBenchmarkGrade g=new RrnaSeedBenchmarkGrade(family,negative);g.loadTruth(truth);g.loadCalls(calls);g.grade();g.write(prefix,arm,clade,assembly);
	}

	public RrnaSeedBenchmarkGrade(String f,boolean neg){family=f;negative=neg;}
	public Bucket bucket(String seq){
		assert(seq!=null&&!seq.isEmpty()):"Every grading bucket is one contig in one explicit assembly";
		Bucket b=bySeq.get(seq);if(b==null){b=new Bucket(seq);bySeq.put(seq,b);ordered.add(b);}return b;
	}
	public void loadTruth(String path){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> keys=new ObjectSet<String>(String.class);
		try{final byte[] h=in.nextLine();require(h!=null&&new String(h,StandardCharsets.UTF_8).equals("seqid\tfamily\tstart\tstop\tstrand"),"Expected five-column CM truth schema");
			for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==5&&p.termEquals("5S",1),"Only the reviewed 5S CM truth is accepted");
				final Locus l=new Locus(p.parseString(0),family,p.parseInt(2),p.parseInt(3),strand(p.parseString(4)));
				require(keys.add(l.seqid+"\t"+l.start+"\t"+l.stop+"\t"+l.strand),"Duplicate CM truth coordinates");bucket(l.seqid).truth.add(l);referenceCount++;}
		}finally{require(!in.close(),"CM truth read failed");}
	}
	public void loadCalls(String path){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){if(row.length==0||row[0]=='#'){continue;}p.set(row);require(p.terms()==9,"Expected nine GFF columns");
			if(!p.termEquals("rRNA",2)&&!p.termEquals("RNA",2)){continue;}final String f=family(row,p);
			require(f.equals("18S")||f.equals("euk5S")||f.equals("16S")||f.equals("23S")||f.equals("prok5S"),"Unexpected enabled RNA family: "+f);
			require(family.equals("euk5S")?(f.equals("18S")||f.equals("euk5S")):(!f.equals("18S")&&!f.equals("euk5S")),"Mixed domain profile in integrated calls: "+f);
			rnaRows++;if(!f.equals(family)){otherFamilyRows++;continue;}
			final Locus l=new Locus(p.parseString(0),f,p.parseInt(3),p.parseInt(4),strand(p.parseString(6)));bucket(l.seqid).calls.add(l);callCount++;}
		}finally{require(!in.close(),"Integrated GFF read failed");}
		require(rnaRows==callCount+otherFamilyRows,"RNA family partition must conserve all caller rows");
	}
	public static String family(byte[] row,LineParser1 p){
		return family(row,p,true);
	}
	/** Legacy CallGenes rRNA has a type token but no family tag; only explicit mixed-mode readers allow it. */
	static String family(byte[] row,LineParser1 p,boolean required){
		assert(p.terms()==9):"Family attribution reads the GFF attribute column";
		p.setBounds(8);final int a=p.a(),end=p.b();String result=null;final String tag="family:";
		for(int pos=a;pos<end;pos++){
			if(pos!=a && row[pos-1]!=',' && row[pos-1]!=';'){continue;}if(pos+tag.length()>end){continue;}boolean match=true;
			for(int i=0;i<tag.length();i++){if(row[pos+i]!=(byte)tag.charAt(i)){match=false;break;}}if(!match){continue;}
			int stop=pos+tag.length();while(stop<end && row[stop]!=',' && row[stop]!=';' && !Character.isWhitespace((char)row[stop])){stop++;}
			require(result==null&&stop>pos+tag.length(),"Missing/duplicate family value");result=new String(row,pos+tag.length(),stop-pos-tag.length(),StandardCharsets.US_ASCII);
		}
		require(!required||result!=null,"Integrated RNA call lacks explicit family attribution");return result;
	}
	public void grade(){
		assert(tp==0&&fp==0&&fn==0):"An input pair is graded only once";
		for(Bucket b:ordered){
			b.match=NcrnaGffFamilyGrader.matchOneToOne(b.calls,negative?new ArrayList<Locus>():b.truth,0.5f);
			tp+=b.match.matched;fp+=b.calls.size()-b.match.matched;if(!negative){fn+=b.truth.size()-b.match.matched;}
			final LongHashSet plus=new LongHashSet(),minus=new LongHashSet();
			for(int i=0;i<b.calls.size();i++){
				final Locus call=b.calls.get(i);final long key=((long)call.start<<32)|(call.stop&0xffffffffL);
				// Exact query repetition is independent of truth overlap or the matcher outcome.
				if(!(call.strand=='+'?plus:minus).add(key)){duplicates++;}final int j=b.match.callMatchedTruthIdx[i];
				if(j>=0){final MatchedPair pair=new MatchedPair(family,b.truth.get(j),b.calls.get(i));if(pair.d5==0){exact5++;}if(pair.d3==0){exact3++;}if(pair.d5==0&&pair.d3==0){both++;}}
				if(negative&&overlapAny(b.calls.get(i),b.truth)){controlOverlap++;}
			}
		}
		require(tp+fp==callCount&&(negative?tp==0&&fn==0:tp+fn==referenceCount),"One-to-one grading must conserve calls and biological truth");
	}
	static boolean overlapAny(Locus call,ArrayList<Locus> truth){
		assert(call!=null&&truth!=null):"Control overlap is a spatial diagnostic, never a true positive";
		for(Locus t:truth){if(call.seqid.equals(t.seqid)&&Math.min(call.stop,t.stop)>=Math.max(call.start,t.start)){return true;}}return false;
	}
	public void write(String prefix,String arm,String clade,String assembly){
		assert(tp+fp==callCount):"Write only a completed grade";
		RrnaCmWindowParity.fresh(prefix,new String[]{".summary.tsv",".loci.tsv.gz"});
		final ByteStreamWriter summary=RrnaSourceTargetPrep.writer(prefix+".summary.tsv"),loci=RrnaSourceTargetPrep.writer(prefix+".loci.tsv.gz");final ByteBuilder b=new ByteBuilder();
		try{summary.println("arm\tclade\tassembly\tfamily\tnegative\ttruth\tcontrolTruth\tcalls\tTP\tFP\tFN\tprecision\trecall\texact5\texact3\tbothExact\texact5Fraction\texact3Fraction\tbothExactFraction\tduplicateCalls\tcontrolOverlapAnyStrand\totherFamilyCalls");
			b.append(arm).tab().append(clade).tab().append(assembly).tab().append(family).tab().append(negative?1:0).tab().append(negative?0:referenceCount).tab().append(negative?referenceCount:0).tab().append(callCount).tab().append(tp).tab().append(fp).tab().append(fn).tab();
			ratio(b,tp,tp+fp);b.tab();ratio(b,tp,tp+fn);b.tab().append(exact5).tab().append(exact3).tab().append(both).tab();ratio(b,exact5,tp);b.tab();ratio(b,exact3,tp);b.tab();ratio(b,both,tp);b.tab().append(duplicates).tab().append(controlOverlap).tab().append(otherFamilyRows).nl();summary.print(b);
			loci.println("assembly\tseqid\trole\tindex\tstart\tstop\tstrand\toutcome\tpartnerIndex\td5\td3\tcontrolOverlapAnyStrand");
			for(Bucket bucket:ordered){for(int i=0;i<bucket.calls.size();i++){writeLocus(loci,b,assembly,bucket,true,i);}for(int i=0;i<bucket.truth.size();i++){writeLocus(loci,b,assembly,bucket,false,i);}}
		}finally{RrnaCmWindowParity.finish(summary,loci);}
	}
	void writeLocus(ByteStreamWriter out,ByteBuilder b,String assembly,Bucket bucket,boolean query,int i){
		assert(bucket.match!=null):"Locus partners come from the completed maximum-cardinality matching";
		final Locus l=(query?bucket.calls:bucket.truth).get(i);final int partner=negative?-1:(query?bucket.match.callMatchedTruthIdx[i]:bucket.match.truthMatchedCallIdx[i]);
		b.clear().append(assembly).tab().append(l.seqid).tab().append(query?"CALL":negative?"CONTROL_TRUTH":"TRUTH").tab().append(i).tab().append(l.start).tab().append(l.stop).tab().append(l.strand).tab().append(partner>=0?"TP":query?"FP":negative?"CONTROL":"FN").tab().append(partner).tab();
		if(partner>=0){final MatchedPair p=new MatchedPair(family,query?bucket.truth.get(partner):l,query?l:bucket.calls.get(partner));b.append(p.d5).tab().append(p.d3);}else{b.append("NA\tNA");}
		b.tab().append(query&&negative&&overlapAny(l,bucket.truth)?1:0).nl();out.print(b);
	}
	static void ratio(ByteBuilder b,long numerator,long denominator){assert(numerator>=0&&denominator>=numerator):"A fraction counts a subset of its explicit denominator";if(denominator==0){b.append("NA");}else{b.append(numerator/(double)denominator,8);}}
	public static char strand(String s){require(s.equals("+")||s.equals("-"),"Known strand required");return s.charAt(0);}
	static void require(boolean b,String reason){if(!b){throw new IllegalArgumentException(reason);}}
	public static final class Bucket{Bucket(String s){seq=s;}final String seq;public final ArrayList<Locus> calls=new ArrayList<Locus>(),truth=new ArrayList<Locus>();MatchResult match;}
	final String family;public final boolean negative;
	final ObjectMap<String,Bucket> bySeq=new ObjectMap<String,Bucket>(String.class,Bucket.class);final ArrayList<Bucket> ordered=new ArrayList<Bucket>();
	public long referenceCount,callCount,rnaRows,otherFamilyRows,tp,fp,fn,duplicates,exact5,exact3,both,controlOverlap;
}
