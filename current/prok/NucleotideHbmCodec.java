package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;

/** Lossless compact storage for the nucleotide HBM1 five-count projection.
 * Preserves model headers (including depth), original sequences, model order and
 * every A/C/G/T/deletion count. It is not the amino-acid AAGraph HBMC contract.
 * Rows store A48 counts in pivot-first order, then remaining ACGT, then deletion;
 * trailing zero fields are omitted. A footer makes truncation between models fail.
 * @author Keqing
 */
public final class NucleotideHbmCodec {
	public static void main(String[] args){
		String in=null, out=null;
		for(String arg:Parser.parseConfig(args)){
			final int eq=arg.indexOf('=');
			if(eq<1){throw bad("Expected in=source out=fresh.hbmc[.gz] or out=fresh.hbm[.gz]");}
			final String key=arg.substring(0,eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in") && in==null){in=value;}
			else if(key.equals("out") && out==null){out=value;}
			else{throw bad("Unknown or duplicate option: "+key);}
		}
		if(in==null || out==null || !Tools.testOutputFiles(false,false,false,out)){throw bad("Require input and a fresh output");}
		convert(in,out);
	}

	/** Converts canonical HBM1 text or nucleotide HBMC, without reconstructing metadata. */
	public static void convert(String in,String out){
		if(in==null || out==null){throw bad("Conversion requires input and output paths");}
		final long models,positions;
		try(Reader reader=new Reader(in); Writer writer=new Writer(out,isCompact(out),false)){
			for(Model model=reader.next();model!=null;model=reader.next()){writer.add(model);}
			writer.finish();
			models=writer.models;positions=writer.positions;
		}
		System.err.println("NUCLEOTIDE_HBM_CONVERT_PASS models="+models+" positions="+positions);
	}

	/** Loads the same scoring state as TrnaConsensusBuilder's legacy text reader. */
	public static BaseGraph[] loadModels(String path){
		final ArrayList<BaseGraph> graphs=new ArrayList<BaseGraph>();
		try(Reader reader=new Reader(path)){
			for(Model model=reader.next();model!=null;model=reader.next()){
				int end=1;while(end<model.header.length && model.header[end]!='\t'){end++;}
				final String name=new String(model.header,1,end-1,StandardCharsets.UTF_8);
				final byte[] bases=model.bases.clone();Tools.toUpperCase(bases);
				final BaseGraph graph=new BaseGraph(name,bases,null,graphs.size(),0);
				for(int pos=0, offset=0;pos<bases.length;pos++,offset+=5){
					final int[] c=model.counts;
					BaseGraphHelper.setCounts(graph,pos,c[offset],c[offset+1],c[offset+2],c[offset+3],c[offset+4]);
				}
				BaseGraphHelper.initForScoring(graph);graphs.add(graph);
			}
		}
		return graphs.toArray(new BaseGraph[graphs.size()]);
	}

	/** Writes the same five-count graph projection as the existing text writer. */
	public static void writeModels(ArrayList<BaseGraph> graphs,String out,boolean overwrite){
		if(graphs==null){throw bad("A model list is required");}
		try(Writer writer=new Writer(out,true,overwrite)){
			for(BaseGraph graph:graphs){
				if(graph==null || graph.original.length==0){throw bad("Cannot write an empty model");}
				final int[] counts=new int[Math.multiplyExact(graph.original.length,5)];
				for(int i=0, offset=0;i<graph.original.length;i++,offset+=5){
					System.arraycopy(graph.ref[i].acgtCount,0,counts,offset,4);counts[offset+4]=graph.del[i].countSum;
					validateCounts(counts,offset);
				}
				final ByteBuilder header=new ByteBuilder();
				header.append('>').append(graph.name).tab().append("depth=")
					.append((long)graph.ref[0].countSum+graph.del[0].countSum);
				writer.add(new Model(header.toBytes(),graph.original,counts));
			}
			writer.finish();
		}
	}

	public static boolean isCompact(String path){
		if(path==null){return false;}
		final String lower=path.toLowerCase(Locale.ROOT);
		return lower.endsWith(".hbmc") || lower.endsWith(".hbmc.gz");
	}

	/** One model is retained during conversion; graph-loading deliberately retains all graphs. */
	static final class Model {
		Model(byte[] header_,byte[] bases_,int[] counts_){
			if(header_==null || header_.length<2 || header_[0]!='>' || bases_==null || bases_.length<1 ||
					counts_==null || (long)bases_.length*5!=counts_.length){throw bad("Invalid model shape");}
			header=header_;bases=bases_;counts=counts_;
		}
		final byte[] header,bases;
		final int[] counts;
	}

	/** Strict framing rejects short rows, early model boundaries and missing compact footers. */
	static final class Reader implements AutoCloseable {
		Reader(String path){
			input=ByteFile.makeByteFile(path,false);
			final byte[] first=input.nextLine();
			compact=Arrays.equals(first,COMPACT);
			if(!compact && !Arrays.equals(first,DENSE)){input.close();throw bad("Unrecognized nucleotide HBM header: "+path);}
		}
		Model next(){
			if(ended){return null;}
			final byte[] header=input.nextLine();
			if(header==null){
				if(compact){throw bad("Missing nucleotide HBMC completion footer");}
				ended=true;return null;
			}
			if(compact && Tools.startsWith(header,"#end\t")){
				final String expected="#end\t"+models+"\t"+positions;
				if(!Arrays.equals(header,expected.getBytes(StandardCharsets.US_ASCII)) || input.nextLine()!=null){throw bad("Invalid completion footer or trailing data");}
				ended=true;return null;
			}
			if(header.length<2 || header[0]!='>'){throw bad("Expected model header at model "+models);}
			final byte[] bases=input.nextLine();
			if(bases==null || bases.length<1 || bases.length>MAX_LENGTH){throw bad("Missing or excessive model sequence at model "+models);}
			for(byte b:bases){if(b<33 || b>126){throw bad("Sequence must contain printable ASCII without whitespace");}}
			final int[] counts=new int[bases.length*5];
			for(int pos=0;pos<bases.length;pos++){
				final byte[] row=input.nextLine();
				if(row==null){throw bad("Truncated counts at model "+models+", position "+pos);}
				decode(row,compact,pivot(bases[pos]),counts,pos*5);
			}
			models++;positions+=bases.length;
			return new Model(header,bases,counts);
		}
		@Override public void close(){if(input.close()){throw bad("Nucleotide HBM input close failed");}}
		final ByteFile input;
		final boolean compact;
		boolean ended;
		long models,positions;
	}

	/** Uses the standard BBTools writer; filename controls the outer gzip container. */
	static final class Writer implements AutoCloseable {
		Writer(String path,boolean compact_,boolean overwrite){
			compact=compact_;output=new ByteStreamWriter(path,overwrite,false,false);output.start();
			buffer.append(compact?COMPACT:DENSE).nl();flush();
		}
		void add(Model model){
			assert(!finished) : "No model may be appended after the completion footer";
			buffer.append(model.header).nl().append(model.bases).nl();flush();
			for(int pos=0;pos<model.bases.length;pos++){
				final int offset=pos*5;validateCounts(model.counts,offset);
				if(compact){encode(model.counts,offset,pivot(model.bases[pos]),buffer);}
				else{for(int c=0;c<5;c++){if(c>0){buffer.tab();}buffer.append(model.counts[offset+c]);}}
				buffer.nl();if(buffer.length()>16384){flush();}
			}
			flush();models++;positions+=model.bases.length;
		}
		void finish(){
			if(finished){throw bad("Writer already finished");}
			if(compact){buffer.append("#end\t").append(models).tab().append(positions).nl();flush();}
			finished=true;
		}
		void flush(){output.print(buffer);buffer.clear();}
		@Override public void close(){if(output.poisonAndWait()){throw bad("Nucleotide HBM output close failed");}}
		final ByteStreamWriter output;
		final boolean compact;
		final ByteBuilder buffer=new ByteBuilder(32768);
		boolean finished;
		long models,positions;
	}

	/** Counts use the canonical A48 primitive after explicit digit and overflow checks. */
	static void decode(byte[] row,boolean compact,int pivot,int[] counts,int offset){
		assert(counts!=null && offset>=0 && offset+5<=counts.length) : "Five count slots are required for each reference position";
		int field=0,start=0;
		for(int end=0;end<=row.length;end++){
			if(end<row.length && row[end]!='\t'){continue;}
			if(field>=5 || end==start){throw bad("Expected one to five nonempty compact counts, or five dense counts");}
			long value=0;
			for(int i=start;i<end;i++){
				final int digit=row[i]-(byte)'0';
				if(digit<0 || digit>=(compact?64:10)){throw bad("Invalid count digit");}
				value=value*(compact?64:10)+digit;
				if(value>Integer.MAX_VALUE){throw bad("Count exceeds signed 32-bit range");}
			}
			final int count=(int)(compact?Parse.parseLongA48(row,start,end):value);
			counts[offset+(compact?column(field,pivot):field)]=count;
			field++;start=end+1;
		}
		if(!compact && field!=5){throw bad("Dense HBM rows require exactly five counts");}
		validateCounts(counts,offset);
	}
	static void encode(int[] counts,int offset,int pivot,ByteBuilder out){
		assert(counts!=null && offset>=0 && offset+5<=counts.length) : "Encoding requires all five positional counts";
		int last=4;while(last>0 && counts[offset+column(last,pivot)]==0){last--;}
		for(int c=0;c<=last;c++){if(c>0){out.tab();}out.appendA48(counts[offset+column(c,pivot)]);}
	}
	private static int column(int field,int pivot){
		assert(field>=0 && field<5 && pivot>=0 && pivot<4) : "Invalid compact count column";
		return field==4?4:field==0?pivot:field<=pivot?field-1:field;
	}
	private static int pivot(byte base){
		final int x=AminoAcid.baseToNumber[base];return x>=0 && x<4?x:0;
	}
	private static void validateCounts(int[] counts,int offset){
		long sum=0;
		for(int c=0;c<5;c++){final int n=counts[offset+c];if(n<0){throw bad("Negative nucleotide count");}sum+=n;}
		if(sum>Integer.MAX_VALUE){throw bad("Positional depth exceeds BaseGraph signed 32-bit sums");}
	}
	private static IllegalArgumentException bad(String why){return new IllegalArgumentException(why);}
	private static final byte[] DENSE="#HBM\t1".getBytes(StandardCharsets.US_ASCII);
	private static final byte[] COMPACT="#HBMC_NT\t1".getBytes(StandardCharsets.US_ASCII);
	private static final int MAX_LENGTH=1000000;
	private NucleotideHbmCodec(){}
}
