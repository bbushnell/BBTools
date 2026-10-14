package assemble;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Exact-byte sampling fixtures with interleaved held-out provenance. @author Fischl */
public final class SameKTrainingSampleTest {

	/** Generates a source fixture or checks an independently emitted subset. */
	public static void main(final String[] args){
		if(args.length==2 && args[0].equals("generate")){generate(args[1]);}
		else if(args.length==4 && args[0].equals("check")){check(args[1], args[2], Integer.parseInt(args[3]));}
		else{throw new IllegalArgumentException("generate NEW_directory, or check input_directory output_directory quota");}
	}

	/** Writes twelve candidates, forty controls, and irrelevant held-out metadata between them. */
	private static void generate(final String root){
		if(!new File(root).mkdir()){throw new IllegalArgumentException("Fixture directory must be new.");}
		final ByteStreamWriter data=new ByteStreamWriter(root+"/train.tsv", false, false, false);
		final ByteStreamWriter rows=new ByteStreamWriter(root+"/rows.tsv", false, false, false);
		data.start(); rows.start();
		data.print("#dims\t27\t1\n#schema=samek27-v1; unknown labels excluded; no class balancing\n");
		rows.print("fold\trow\tgenome\tkind\tevent\ttarget\n");
		final ByteBuilder vector=new ByteBuilder(), meta=new ByteBuilder();
		for(int i=1; i<=52; i++){
			vector.clear();
			for(int j=0; j<27; j++){if(j>0){vector.tab();} vector.append(i+j, 3);}
			data.print(vector.tab().append(i%2).nl());
			rows.print(meta.clear().append("train\t").append(i).append("\tg0001\t")
					.append(i<=12 ? "candidate" : "control").tab().append(i<=12 ? i : i-12).tab().append(i%2).nl());
			rows.print(meta.clear().append("test\t").append(i).append("\tg0002\tcandidate\t").append(i).append("\t1\n"));
		}
		if(data.poisonAndWait() | rows.poisonAndWait()){throw new IllegalStateException("Fixture write failed.");}
	}

	/** Checks retained source bytes, target identity, stable order and complete candidate retention. */
	private static void check(final String input, final String output, final int quota){
		final ArrayList<byte[]> source=new ArrayList<byte[]>();
		final ByteFile original=ByteFile.makeByteFile(input+"/train.tsv", true);
		try{
			original.nextLine(); original.nextLine();
			for(byte[] line=original.nextLine(); line!=null; line=original.nextLine()){source.add(line);}
		}finally{if(original.close()){throw new AssertionError("Cannot read source fixture.");}}
		final ByteFile data=ByteFile.makeByteFile(output+"/train.tsv", true), meta=ByteFile.makeByteFile(output+"/rows.tsv", true);
		final LineParser1 lp=new LineParser1('\t');
		int seen=0, candidates=0, controls=0, previous=0;
		try{
			data.nextLine(); data.nextLine(); meta.nextLine();
			for(byte[] line=meta.nextLine(); line!=null; line=meta.nextLine()){
				lp.set(line);
				final int index=lp.parseInt(1);
				if(lp.terms()!=6 || lp.parseInt(0)!=++seen || index<=previous || index>52 ||
						!lp.termEquals("g0001", 2) || lp.parseInt(5)!=index%2 || !Arrays.equals(data.nextLine(), source.get(index-1))){
					throw new AssertionError("Subset bytes/identity/order differ from source fixture.");
				}
				previous=index;
				if(index<=12){
					if(!lp.termEquals("candidate", 3) || index!=++candidates){throw new AssertionError("Candidate omitted or mislabeled.");}
				}else{
					if(!lp.termEquals("control", 3)){throw new AssertionError("Control mislabeled.");}
					controls++;
				}
			}
			if(candidates!=12 || controls!=quota || data.nextLine()!=null){throw new AssertionError("Subset count/EOF differs.");}
		}finally{if(data.close() | meta.close()){throw new AssertionError("Cannot read subset fixture.");}}
		System.err.println("SAME_K_SAMPLE_FIXTURE_PASS candidates="+candidates+" controls="+controls);
	}
}
