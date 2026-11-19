package ddl;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Objects;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Exact K5/K6 occurrence census for one explicitly declared protein population.
 * Real and shadow populations are counted separately. This producer does not
 * choose their mixture, inspect evaluation data, or enable a blacklist in a model.
 * Input identities are inherited labels, never recomputed digests.
 * @author Yaoyao
 */
public final class ProteinKmerFrequency{

	public static void main(String[] args) throws Exception{
		final ProteinKmerFrequency tool=new ProteinKmerFrequency(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	/** Standard argument parsing with an exact record census and explicit length domain. */
	public ProteinKmerFrequency(String[] args){
		final PreParser pp=new PreParser(args,getClass(),false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		final HashSet<String> seen=new HashSet<String>();
		for(String arg:pp.args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			final String key=pair[0].toLowerCase(Locale.ROOT),value=pair.length>1 ? pair[1] : null;
			require(seen.add(key),"Duplicate option: "+key);
			if(key.equals("k")){k=Integer.parseInt(value);}
			else if(key.equals("expected")){expected=Parse.parseKMG(value);}
			else if(key.equals("minaa")){minAA=Parse.parseIntKMG(value);}
			else if(key.equals("maxaa")){maxAA=Parse.parseIntKMG(value);}
			else if(key.equals("top")){top=Parse.parseIntKMG(value);}
			else if(key.equals("population")){population=value;}
			else if(key.equals("insha80")){inheritedPin=value;}
			else if(key.equals("t") || key.equals("threads")){
				require(Integer.parseInt(value)==1,"The dense counter uses one worker; use t=1");
			}else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown option: "+arg);}
		}
		require(k==5 || k==6,"Require k=5 or k=6");
		require(expected>0 && minAA>0 && maxAA>=minAA,"Require expected>0 and explicit minaa/maxaa");
		require(top>0 && top<=universe(k),"top must fit the word universe");
		require("real".equals(population) || "shadow".equals(population),"Require population=real or shadow");
		require(inheritedPin!=null && (inheritedPin.matches("[0-9a-f]{20}") || inheritedPin.equals("unmeasured")),
			"Require inherited insha80 or explicit unmeasured; this producer computes no hashes");
		require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null &&
			!parser.append && !parser.overwrite && parser.maxReads<0,"Require one complete input and a new output directory");
		input=Paths.get(parser.in1);output=Paths.get(parser.out1);
	}

	/** Stream once, conserve every valid window, and publish PASS only after all writers close. */
	private void process() throws Exception{
		require(!Files.exists(output,LinkOption.NOFOLLOW_LINKS),"Output already exists: "+output);
		final BasicFileAttributes before=attributes(input);
		Shared.setThreads(1);Shared.AMINO_IN=true;
		ReadWrite.setZipThreads(1);
		ReadWrite.ZIPLEVEL=9;ReadWrite.FORCE_BGZIP=ReadWrite.USE_BGZIP=true;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
		final FileFormat format=FileFormat.testInput(input.toString(),FileFormat.FASTA,null,true,false);
		require(format.fasta() && format.amino(),"Require protein FASTA (.faa or amino input mode)");
		final long[] counts=new long[universe(k)];
		final Streamer reader=StreamerFactory.makeStreamer(format,null,0,true,-1,false,true,1);
		long records=0,residues=0,windows=0;
		reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList();batch!=null && batch.size()>0;batch=reader.nextList()){
				for(Read read:batch){
					require(read.mate==null && read.id!=null && !read.id.isEmpty() && read.bases!=null,
						"Malformed or paired protein record");
					require(read.bases.length>=minAA && read.bases.length<=maxAA,"Protein outside declared length domain: "+read.id);
					records=Math.addExact(records,1);
					require(records<=expected,"Input has more records than expected");
					residues=Math.addExact(residues,read.bases.length);
					windows=Math.addExact(windows,count(read.bases,k,counts));
				}
				if(records>=nextProgress){outstream.println("Counted "+records+" records and "+windows+" valid K"+k+" windows.");nextProgress=records+1000000;}
			}
		}finally{reader.close();}
		require(!reader.errorState() && records==expected,"Input read failed or record count differs from expected: "+records+" != "+expected);
		unchanged(before,attributes(input));
		long sum=0,distinct=0;
		for(long n:counts){sum=Math.addExact(sum,n);if(n>0){distinct++;}}
		require(sum==windows && windows>0,"Count conservation failed or no valid words");
		final int[] selected=select(counts,top);
		Files.createDirectory(output);
		writeCounts(counts);
		writeTop(counts,selected);
		final String summary="population\tk\trecords\tresidues\tvalid_words\tdistinct_words\tselected_words\n"+
			population+"\t"+k+"\t"+records+"\t"+residues+"\t"+windows+"\t"+distinct+"\t"+selected.length+"\n";
		writeText("summary.tsv",summary);
		writeText("recipe.txt","in="+input+"\ninsha80="+inheritedPin+"\nidentity_verification=inherited_not_recomputed\n"+
			"bulk_hashing=false\npopulation="+population+"\nk="+k+"\nexpected="+expected+"\nminaa="+minAA+"\nmaxaa="+maxAA+
			"\nalphabet="+AA+"\ncount=occurrences\nambiguous_residue=break_window\nrecord_boundary=break_window\n"+
			"trimming=none\ntop="+top+"\nties=ascending_packed_word\nnew_output_identity=unmeasured\ncompression=BGZF9\n");
		unchanged(before,attributes(input));
		writeText("PASS","PROTEIN_KMER_FREQUENCY_PASS\n");
		outstream.println("PROTEIN_KMER_FREQUENCY_PASS population="+population+" k="+k+" records="+records+" windows="+windows);
	}

	/** Count literal forward words. Case is ignored; every nonstandard residue resets the window. */
	public static long count(byte[] peptide,int k,long[] counts){
		require(peptide!=null && counts!=null && counts.length==universe(k),"Invalid peptide or dense count universe");
		final int modulus=counts.length/20;
		int word=0,valid=0;long windows=0;
		for(byte b:peptide){
			final int aa=b>=0 && b<ENCODE.length ? ENCODE[b] : -1;
			if(aa<0){word=valid=0;continue;}
			word=(word%modulus)*20+aa;
			if(valid<k){valid++;}
			if(valid==k){counts[word]=Math.addExact(counts[word],1);windows++;}
		}
		return windows;
	}

	/** Keep positive-count words only. A bounded worst-first heap avoids sorting the whole universe. */
	public static int[] select(long[] counts,int limit){
		require(counts!=null && limit>0 && limit<=counts.length,"Invalid top-list request");
		final int[] heap=new int[limit];int size=0;
		for(int word=0;word<counts.length;word++){
			require(counts[word]>=0,"Negative occurrence count");
			if(counts[word]==0){continue;}
			if(size<limit){
				int child=size++;heap[child]=word;
				while(child>0){
					final int parent=(child-1)/2;
					if(!worse(heap[child],heap[parent],counts)){break;}
					final int swap=heap[parent];heap[parent]=heap[child];heap[child]=swap;child=parent;
				}
			}else if(worse(heap[0],word,counts)){heap[0]=word;sink(heap,size,counts);}
		}
		final int[] result=new int[size];
		for(int i=size-1;i>=0;i--){result[i]=heap[0];heap[0]=heap[--size];if(size>0){sink(heap,size,counts);}}
		return result;
	}

	private static boolean worse(int a,int b,long[] counts){return counts[a]<counts[b] || (counts[a]==counts[b] && a>b);}
	private static void sink(int[] heap,int size,long[] counts){
		int parent=0;
		for(int child=1;child<size;child=parent*2+1){
			if(child+1<size && worse(heap[child+1],heap[child],counts)){child++;}
			if(!worse(heap[child],heap[parent],counts)){break;}
			final int swap=heap[parent];heap[parent]=heap[child];heap[child]=swap;parent=child;
		}
	}

	/** Dense base20 indexes have the same ordering as the canonical five-bit packed words. */
	public static int packed(int index,int k){
		require(index>=0 && index<universe(k),"Word index outside declared universe");
		int packed=0;
		for(int shift=0;shift<5*k;shift+=5){packed|=(index%20)<<shift;index/=20;}
		return packed;
	}
	public static String word(int index,int k){
		final int packed=packed(index,k);final char[] letters=new char[k];
		for(int i=0;i<k;i++){letters[i]=AA.charAt((packed>>>(5*(k-i-1)))&31);}
		return new String(letters);
	}
	public static int universe(int k){
		require(k==5 || k==6,"Require k=5 or k=6");return k==5 ? 3200000 : 64000000;
	}

	private void writeCounts(long[] counts){
		final ByteStreamWriter writer=writer("counts.tsv.gz");final ByteBuilder row=new ByteBuilder();
		writer.println("packed_word\tword\toccurrences");
		try{
			for(int i=0;i<counts.length;i++){
				if(counts[i]>0){
					final int packed=packed(i,k);row.clear().append(packed).tab();
					for(int shift=5*(k-1);shift>=0;shift-=5){row.append(AA.charAt((packed>>>shift)&31));}
					writer.print(row.tab().append(counts[i]).nl());
				}
			}
		}finally{require(!writer.poisonAndWait(),"Count table write failed");}
	}
	private void writeTop(long[] counts,int[] selected){
		final ByteStreamWriter fasta=writer("top.faa.gz"),table=writer("top.tsv.gz");
		final ByteBuilder row=new ByteBuilder();table.println("rank\tpacked_word\tword\toccurrences");
		try{
			for(int i=0;i<selected.length;i++){
				final int index=selected[i];final String word=word(index,k);
				fasta.print(row.clear().append('>').append(i+1).append("_count_").append(counts[index]).nl().append(word).nl());
				table.print(row.clear().append(i+1).tab().append(packed(index,k)).tab().append(word).tab().append(counts[index]).nl());
			}
		}finally{require(!(fasta.poisonAndWait() | table.poisonAndWait()),"Top-list write failed");}
	}
	private ByteStreamWriter writer(String name){
		final ByteStreamWriter writer=new ByteStreamWriter(FileFormat.testOutput(output.resolve(name).toString(),FileFormat.TXT,null,false,false,false,false));
		writer.start();return writer;
	}
	private void writeText(String name,String text) throws Exception{Files.write(output.resolve(name),text.getBytes(StandardCharsets.UTF_8));}
	private static BasicFileAttributes attributes(Path path) throws Exception{
		final BasicFileAttributes a=Files.readAttributes(path,BasicFileAttributes.class,LinkOption.NOFOLLOW_LINKS);
		require(a.isRegularFile(),"Input must be a regular nonsymlink file: "+path);return a;
	}
	private static void unchanged(BasicFileAttributes a,BasicFileAttributes b){
		require(a.size()==b.size() && a.lastModifiedTime().equals(b.lastModifiedTime()) && Objects.equals(a.fileKey(),b.fileKey()),"Input changed during census");
	}
	private static int[] encoding(){
		final int[] codes=new int[128];Arrays.fill(codes,-1);
		for(int i=0;i<AA.length();i++){codes[AA.charAt(i)]=codes[Character.toLowerCase(AA.charAt(i))]=i;}
		return codes;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}

	private final PrintStream outstream;
	private Path input,output;
	private String population,inheritedPin;
	private int k=-1,minAA=-1,maxAA=-1,top=2048;
	private long expected=-1,nextProgress=1000000;
	private static final String AA="ARNDCEQGHILKMFPSTWYV";
	private static final int[] ENCODE=encoding();
}
