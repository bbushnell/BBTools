package synth;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;

/** Independent string-graph oracle and exact output-partition audit for branch shredding.
 * String allocation is confined to tiny fixtures; real-genome audit compares the
 * actual source and piece bytes and reports an ordinary length distribution.
 * @author Brian Bushnell, Nilou */
public final class ShredBranchTest {

	/** Runs tiny fixtures, or audits saved shreds without changing their bytes. */
	public static void main(String[] args) throws Exception{
		String mode="fixture", in=null, shreds=null, out=null, name="assembly";
		int k=0;
		final Parser parser=new Parser();
		for(String arg:args){
			final int split=arg.indexOf('=');
			final String key=(split<0 ? arg : arg.substring(0, split)).toLowerCase(Locale.ROOT);
			final String value=split<0 ? null : arg.substring(split+1);
			if(key.equals("mode")){mode=value;}
			else if(key.equals("in")){in=value;}
			else if(key.equals("shreds")){shreds=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("name")){name=value;}
			else if(key.equals("k")){k=Parse.parseIntKMG(value);}
			else if(!parser.parse(arg,key,value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		Parser.processQuality();
		if(out==null){throw new IllegalArgumentException("Require out=fresh directory for fixtures or out=fresh TSV for audit");}
		if(mode.equals("fixture")){
			Shared.setThreads(1);
			final Path root=Paths.get(out); Files.createDirectory(root);
			fixture(root, "repeat", new String[]{"AAAAAA","AAAAAA"}, 3, false, -1);
			final ArrayList<Read> fork=fixture(root, "forward", new String[]{"AAAC","AAAG"}, 3, false, -1);
			require(fork.size()==4 && new String(fork.get(0).bases,java.nio.charset.StandardCharsets.US_ASCII).equals("AAA")
				&& fork.get(1).length()==1, "Literal forward fork must cut after the branch stop");
			fixture(root, "backward", new String[]{"CAAGCCTA","GAAGCCTA"}, 5, false, -1);
			fixture(root, "reverse", new String[]{"GTTT","CTTT"}, 3, false, -1);
			fixture(root, "ambiguity", new String[]{"NAAACNNNAAAGNN","AC","N"}, 3, true, -1);
			fixture(root, "limit", new String[]{"AAAAAA","AAAAAA","AAAC"}, 3, false, 2);
			final String[] dense=new String[65];
			for(int i=0; i<64; i++){dense[i]=""+BASES[i>>4]+BASES[(i>>2)&3]+BASES[i&3];}
			dense[64]="ACGTACGTACGT";
			final ArrayList<Read> shortPieces=fixture(root, "dense", dense, 3, false, -1);
			int single=0; for(Read read:shortPieces){if(read.length()==1){single++;}}
			require(single>=9, "Consecutive branch k-mers must retain sub-k pieces");
			final Random random=new Random(7421);
			final StringBuilder sequence=new StringBuilder();
			for(int i=0; i<600; i++){sequence.append(BASES[random.nextInt(4)]);}
			final String shared=sequence.toString();
			for(int size:new int[]{1,2,4,21,31,32,33,51,63,64,65,127,255,511}){
				fixture(root, "k"+size, new String[]{"A"+shared+"C","G"+shared+"T",reverseComplement("A"+shared+"C")},size,false,-1);
			}
			for(String option:new String[]{"length=4","minlen=2","maxlen=10","median=5","variance=2","mode=linear","equal=t","overlap=1","increment=1","maxns=0"}){
				boolean rejected=false;
				try{new Shred(new String[]{"in="+root.resolve("repeat.fa"),"k=3",option});}
				catch(IllegalArgumentException expected){rejected=true;}
				require(rejected,"Branch mode accepted conflicting option "+option);
			}
			System.out.println("SHRED_BRANCH_FIXTURES_PASS string_graph_oracle=true exact_bases_qualities=true k1_to511=true hashonly_exact_bytes=true legacy_flags_rejected=true");
		}else if(mode.equals("audit")){
			if(in==null || shreds==null || k<1){throw new IllegalArgumentException("audit requires in=assembly shreds=pieces k=positive out=new.tsv");}
			audit(in, shreds, k, out, name);
		}else if(mode.equals("storage")){
			storage(out);
		}else{throw new IllegalArgumentException("mode must be fixture, audit or storage");}
	}

	/** Measures actual primitive backing arrays; headers, vacant slots and work buffers are not conflated with key size. */
	private static void storage(String output) throws Exception{
		final ByteBuilder text=new ByteBuilder("k\tmode\tentries\tallocated_slots\tkey_array_bytes\tcount_array_bytes\tbytes_per_slot\tallocated_array_bytes_per_entry\n");
		for(int k:new int[]{51,127,255,511}){
			for(boolean hashed:new boolean[]{false,true}){
				final ukmer.AbstractKmerTableU table=hashed ? new ukmer.HashArrayH1D(new int[]{1009})
					: new ukmer.HashArrayU1D(new int[]{1009},ukmer.Kmer.getK(k),k);
				final Random random=new Random(5831);
				for(int i=0;i<256;i++){
					final ukmer.Kmer word=new ukmer.Kmer(k);
					for(int b=0;b<k;b++){word.addRightNumeric(random.nextInt(4));}
					require(table.incrementAndReturnNumCreated(word)==1 && table.getValue(word)==1,"Storage probe requires distinct retained keys");
				}
				final int[] values=(int[])field(table,hashed ? ukmer.HashArrayH1D.class : ukmer.HashArrayU1D.class,"values");
				long keyBytes=0;
				if(hashed){
					keyBytes=8L*((long[])field(table,ukmer.HashArrayH1D.class,"hash1")).length
						+8L*((long[])field(table,ukmer.HashArrayH1D.class,"hash2")).length;
				}else{
					for(long[] words:(long[][])field(table,ukmer.HashArrayU.class,"arrays")){keyBytes+=8L*words.length;}
				}
				int occupied=0; for(int value:values){if(value>0){occupied++;}}
				require(occupied==256 && table.size()==256,"Probe must fit entirely in main arrays, without victim entries");
				final long countBytes=4L*values.length, slotBytes=(keyBytes+countBytes)/values.length;
				require(slotBytes==(hashed ? 20 : 8L*ukmer.Kmer.getMult(k)+4),"Backing arrays disagree with documented table payload");
				text.append(k).tab().append(hashed ? "hashonly" : "exact").tab().append(occupied).tab().append(values.length)
					.tab().append(keyBytes).tab().append(countBytes).tab().append(slotBytes)
					.tab().appendSlow((keyBytes+countBytes)/(double)occupied).nl();
			}
		}
		write(output,text);
		System.out.println("SHRED_STORAGE_PASS measured_primitive_arrays=true excludes_headers_work_buffers_and_resize=true");
	}

	/** Test-only inspection avoids adding instrumentation to the production table API. */
	private static Object field(Object object, Class<?> owner, String name) throws Exception{
		final java.lang.reflect.Field field=owner.getDeclaredField(name); field.setAccessible(true); return field.get(object);
	}

	/** The oracle enumerates oriented string neighbors in a canonical set, independently of native bit shifts. */
	private static ArrayList<Read> fixture(Path root, String name, String[] sequences, int k, boolean quality, int limit) throws Exception{
		final Path in=root.resolve(name+(quality ? ".fq" : ".fa")), out=root.resolve(name+".out"+(quality ? ".fq" : ".fa"));
		final ByteBuilder text=new ByteBuilder();
		for(int i=0; i<sequences.length; i++){
			text.append(quality ? '@' : '>').append("contig").append(i).nl().append(sequences[i]).nl();
			if(quality){text.append("+\n"); for(int j=0; j<sequences[i].length(); j++){text.append((char)(33+j%41));} text.nl();}
		}
		write(in.toString(),text);
		final String[] args={"in="+in,"out="+out,"k="+k,"reads="+limit,"qin=33","qout=33","t=1"};
		Shred.main(args);
		final Path exact=root.resolve(name+".exact"+(quality ? ".fq" : ".fa"));
		Shred.main(new String[]{"in="+in,"out="+exact,"k="+k,"reads="+limit,"qin=33","qout=33","t=1","hashonly=f"});
		require(Arrays.equals(Files.readAllBytes(out),Files.readAllBytes(exact)), "Default hash-only and full keys differ: "+name);
		final ArrayList<Read> sources=read(in.toString()), actual=read(out.toString());
		final int count=limit<0 ? sources.size() : Math.min(limit,sources.size());
		final HashSet<String> words=new HashSet<String>();
		for(int s=0; s<count; s++){
			final String sequence=sequences[s];
			for(int pos=k; pos<=sequence.length(); pos++){
				final String word=sequence.substring(pos-k,pos);
				if(valid(word)){words.add(canonical(word));}
			}
		}
		int piece=0;
		for(int s=0; s<count; s++){
			int start=0; final String sequence=sequences[s];
			for(int stop=k; stop<sequence.length(); stop++){
				final String word=sequence.substring(stop-k,stop);
				if(valid(word) && branch(word,words)){
					require(piece<actual.size(),"Missing expected branch piece");
					checkPiece(sources.get(s),actual.get(piece++),start,stop); start=stop;
				}
			}
			if(start<sequence.length()){
				require(piece<actual.size(),"Missing tail");
				checkPiece(sources.get(s),actual.get(piece++),start,sequence.length());
			}
		}
		require(piece==actual.size(),"Unexpected pieces after the exact string-graph partition");
		if(name.equals("repeat")){require(actual.size()==2,"Multiplicity alone must not create a branch");}
		return actual;
	}

	/** Checks both orientations without using the production k-mer encoding or neighbor APIs. */
	private static boolean branch(String word, HashSet<String> words){
		int left=0, right=0;
		for(char base:BASES){
			if(words.contains(canonical(base+word.substring(0,word.length()-1)))){left++;}
			if(words.contains(canonical(word.substring(1)+base))){right++;}
		}
		return left>1 || right>1;
	}
	private static boolean valid(String word){
		for(int i=0;i<word.length();i++){if("ACGT".indexOf(word.charAt(i))<0){return false;}}
		return true;
	}
	private static String canonical(String word){final String rc=reverseComplement(word); return word.compareTo(rc)<0 ? word : rc;}
	private static String reverseComplement(String word){
		final char[] bases=new char[word.length()];
		for(int i=0;i<word.length();i++){
			final int code="ACGT".indexOf(word.charAt(i)); require(code>=0,"Oracle requires unambiguous DNA");
			bases[word.length()-i-1]=BASES[3-code];
		}
		return new String(bases);
	}

	/** Verifies the original interval's identity, bases and qualities, not just its total length. */
	private static void checkPiece(Read source, Read piece, int start, int end){
		require(piece.length()==end-start && piece.id.equals(source.id+"_"+start+"-"+(end-1)),"Wrong interval length/header: "+piece.id);
		for(int i=start;i<end;i++){
			require(piece.bases[i-start]==source.bases[i],"Changed original base at "+i);
			if(source.quality!=null){require(piece.quality!=null && piece.quality[i-start]==source.quality[i],"Changed original quality");}
		}
	}

	/** Exact concatenation checks on real inputs; no inference from the shredder's own counters. */
	private static void audit(String in, String shreds, int k, String out, String name){
		final ArrayList<Read> sources=read(in), pieces=read(shreds);
		final IntList lengths=new IntList();
		long inputBp=0, outputBp=0, belowCount=0, belowBp=0;
		int next=0;
		for(Read source:sources){
			inputBp+=source.length();
			for(int start=0; start<source.length();){
				require(next<pieces.size(),"Output ends before the source partition");
				final Read piece=pieces.get(next++); final int end=start+piece.length();
				require(piece.length()>0 && end<=source.length(),"Output overlaps a contig boundary or is empty");
				checkPiece(source,piece,start,end); lengths.add(piece.length()); outputBp+=piece.length();
				if(piece.length()<k){belowCount++; belowBp+=piece.length();} start=end;
			}
		}
		require(next==pieces.size() && inputBp==outputBp && lengths.size>0,"Real output does not partition the input bases exactly");
		Arrays.sort(lengths.array,0,lengths.size);
		long sum=0; int n50=0;
		for(int i=lengths.size-1;i>=0;i--){sum+=lengths.array[i];if(sum>=outputBp/2+outputBp%2){n50=lengths.array[i];break;}}
		final ByteBuilder table=new ByteBuilder("genome\tk\tinput_contigs\tinput_bp\tpieces\toutput_bp\tmin\tp10\tmedian\tp90\tmax\tmean\tn50\tpieces_lt_k\tbp_lt_k\tbp_fraction_lt_k\n");
		table.append(name).tab().append(k).tab().append(sources.size()).tab().append(inputBp).tab().append(lengths.size).tab().append(outputBp)
			.tab().append(lengths.array[0]).tab().append(lengths.array[(lengths.size-1)/10]).tab().append(lengths.array[(lengths.size-1)/2])
			.tab().append(lengths.array[(int)(9L*(lengths.size-1)/10)]).tab().append(lengths.array[lengths.size-1])
			.tab().appendSlow(outputBp/(double)lengths.size).tab().append(n50).tab().append(belowCount).tab().append(belowBp)
			.tab().appendSlow(belowBp/(double)outputBp).nl();
		write(out,table);
		System.out.println("SHRED_PARTITION_AUDIT_PASS bp="+outputBp+" pieces="+lengths.size);
	}
	private static ArrayList<Read> read(String path){
		return StreamerFactory.getReads(-1,false,FileFormat.testInput(path,FileFormat.FASTA,null,true,true),null,null,null);
	}
	private static void write(String path, ByteBuilder text){
		final ByteStreamWriter out=new ByteStreamWriter(path,false,false,true); out.start(); out.print(text);
		if(out.poisonAndWait()){throw new IllegalStateException("Fixture/audit output failed: "+path);}
	}
	private static void require(boolean valid,String message){if(!valid){throw new IllegalStateException(message);}}
	private static final char[] BASES={'A','C','G','T'};
}
