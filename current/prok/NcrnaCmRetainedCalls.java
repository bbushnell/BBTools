package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Exact one-to-one tblout preservation, excluding only the declared E-value
 * rescaling. Every other fixed field, including scores, model spans and duplicate
 * multiplicity, participates. This measures raw rows, not biological loci.
 * @author Ganyu */
public final class NcrnaCmRetainedCalls {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaCmRetainedCalls.class,false);final PrintStream log=pp.outstream;
		String reference=null,candidate=null,model=null,prefix=null;int clen=0;double ga=Double.NaN;final Parser parser=new Parser();
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("reference")){reference=v;}else if(k.equals("candidate")){candidate=v;}else if(k.equals("model")){model=v;}
				else if(k.equals("clen")){clen=Parse.parseIntKMG(v);}else if(k.equals("ga")){ga=Double.parseDouble(v);}else if(k.equals("outprefix")){prefix=v;}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown flag: "+arg);}
			}
			require(reference!=null && candidate!=null && model!=null && prefix!=null && clen>0 && Double.isFinite(ga),"Require reference= candidate= model= clen= ga= outprefix=");
			require(!parser.append && parser.maxReads<0,"Retention requires complete tables and fresh outputs");
			final long[] counts=compare(reference,candidate,model,clen,ga,prefix);
			log.println("CM_RETAINED_CALLS_PASS reference="+counts[0]+" candidate="+counts[1]+" retained="+counts[2]+" reference_only="+counts[3]+" candidate_only="+counts[4]);
		}finally{Shared.closeStream(log);}
	}
	/** Return reference, candidate, matched, reference-only and candidate-only row counts. */
	public static long[] compare(String reference,String candidate,String model,int clen,double ga,String prefix){
		require(prefix!=null && !prefix.isEmpty(),"Fresh output prefix required");RrnaResourceIO.fresh(prefix,".summary.tsv",".rows.tsv");
		final ArrayList<Row> a=load(reference,model,clen,ga),b=load(candidate,model,clen,ga);
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".rows.tsv");final ByteBuilder buffer=new ByteBuilder();
		int i=0,j=0;long matched=0,lost=0,added=0;
		try{
			out.println("status\treference_line\tcandidate_line\treference_evalue\tcandidate_evalue\ttarget\ttarget_accession\tquery\tmodel\tmdl\tmodel_from\tmodel_to\tseq_from\tseq_to\tstrand\ttrunc\tpass\tgc\tbias\tscore\tinc");
			while(i<a.size() || j<b.size()){
				final Row x=i<a.size()?a.get(i):null,y=j<b.size()?b.get(j):null;
				final int cmp=x==null?1:y==null?-1:x.key.compareTo(y.key);
				if(cmp==0){emit(out,buffer,"RETAINED",x,y);i++;j++;matched++;}
				else if(cmp<0){emit(out,buffer,"REFERENCE_ONLY",x,null);i++;lost++;}
				else{emit(out,buffer,"CANDIDATE_ONLY",null,y);j++;added++;}
			}
		}finally{RrnaResourceIO.finish(out);}
		require(matched+lost==a.size() && matched+added==b.size(),"Sorted one-to-one merge must conserve both raw multisets, including duplicates");
		final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		try{summary.println("reference_rows\tcandidate_rows\tretained_except_evalue\treference_only\tcandidate_only\tall_other_fixed_fields_exact");
			summary.print(buffer.clear().append(a.size()).tab().append(b.size()).tab().append(matched).tab().append(lost).tab().append(added).tab().append(lost==0 && added==0).nl());
		}finally{RrnaResourceIO.finish(summary);}
		return new long[]{a.size(),b.size(),matched,lost,added};
	}
	private static ArrayList<Row> load(String path,String model,int clen,double ga){
		require(path!=null && model!=null && clen>0 && Double.isFinite(ga),"Pinned CM identity and metadata required to validate each input table");
		final ArrayList<Row> rows=new ArrayList<Row>();final RrnaCmTruthBuilder.Fields f=new RrnaCmTruthBuilder.Fields();final ByteBuilder key=new ByteBuilder();
		NcrnaMitoCmEvidence.parseTable(path,model,clen,ga,h->{
			f.set(h.raw);key.clear();for(int n=0;n<17;n++){if(n==15){continue;}if(n>0){key.tab();}key.append(h.raw,f.start[n],f.end[n]-f.start[n]);}
			rows.add(new Row(key.toString(),h.evalue,h.line));
		});
		rows.sort(Comparator.comparing(r->r.key));return rows;
	}
	private static void emit(ByteStreamWriter out,ByteBuilder b,String status,Row a,Row z){
		assert(a!=null || z!=null):"Each merge event consumes at least one actual tblout row";
		b.clear().append(status).tab();if(a==null){b.append("NA");}else{b.append(a.line);}b.tab();
		if(z==null){b.append("NA");}else{b.append(z.line);}b.tab().append(a==null?"NA":a.evalue).tab().append(z==null?"NA":z.evalue).tab().append(a==null?z.key:a.key).nl();out.print(b);
	}
	private static final class Row {Row(String k,String e,long n){key=k;evalue=e;line=n;}final String key,evalue;final long line;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
