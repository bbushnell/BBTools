package assemble;

import java.io.File;
import java.nio.charset.StandardCharsets;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Builds controlled model/feature fixtures and checks the native NN trace. @author Fischl */
public final class SameKJoinModelGateTest {

	/** Uses real exported vectors; only isolated negative copies are altered. */
	public static void main(final String[] args){
		if(args.length==3 && args[0].equals("generate")){generate(args[1], args[2]);}
		else if(args.length==3 && args[0].equals("check")){check(args[1], Integer.parseInt(args[2]));}
		else{throw new IllegalArgumentException("generate features.tsv NEW_directory, or check scores.tsv accepted=0|1");}
	}

	/** Creates a constant0.8 boundary model and malformed feature copies without touching their source. */
	private static void generate(final String features, final String root){
		if(!new File(root).mkdir()){throw new IllegalArgumentException("Fixture directory must be new.");}
		for(float bias : new float[]{0.8f, 2}){
			final ByteBuilder model=new ByteBuilder("##bbnet\n#version 1\n#concise\n#dense\n#density 1\n#blocksize 1\n#seed 1\n#layers 2\n#dims 27 1\n#edges 27\n#coding decimal\n##layer 1\nC28 LINEAR ");
			model.appendSlow(bias);
			for(int i=0; i<27; i++){model.append(" 0");}
			write(root+(bias==2 ? "/badscore.bbnet" : "/constant.bbnet"), model.nl());
		}
		for(String variant : new String[]{"badid", "badwidth", "nonfinite", "badk", "missingtrailer"}){
			final ByteFile input=ByteFile.makeByteFile(features, true);
			final ByteStreamWriter out=new ByteStreamWriter(root+"/"+variant+".tsv", false, false, false);
			out.start();
			final LineParser1 fields=new LineParser1('\t');
			final ByteBuilder row=new ByteBuilder();
			boolean changed=false;
			try{
				for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
					final String text=new String(line, StandardCharsets.UTF_8);
					if(variant.equals("missingtrailer") && text.startsWith("#controls")){changed=true; continue;}
					row.clear();
					if(!changed && text.startsWith("candidate\t" ) && !variant.equals("missingtrailer")){
						fields.set(line);
						for(int i=0; i<fields.terms(); i++){
							if(i>0){row.tab();}
							if(i==1 && variant.equals("badid")){row.append(2);}
							else if(i==2 && variant.equals("nonfinite")){row.append("NaN");}
							else if(i==26 && variant.equals("badk")){row.append(0);}
							else{row.append(fields.parseByteArray(i));}
						}
						if(variant.equals("badwidth")){row.append("\t0");}
						changed=true;
					}else{row.append(line);}
					out.print(row.nl());
				}
			}finally{if(input.close() | out.poisonAndWait()){throw new IllegalStateException("Cannot create negative fixture.");}}
			if(!changed){throw new AssertionError("Negative fixture did not mutate its intended field.");}
		}
	}

	/** Checks that only the226 real candidates, not25 virtual controls, entered inference. */
	private static void check(final String path, final int accepted){
		if(accepted!=0 && accepted!=1){throw new IllegalArgumentException("Expected binary decision.");}
		final ByteFile input=ByteFile.makeByteFile(path, true);
		final LineParser1 fields=new LineParser1('\t');
		int event=0;
		try{
			final byte[] header=input.nextLine();
			if(header==null || !new String(header, StandardCharsets.UTF_8).equals("event\tk\tscore_bits\tscore\taccepted")){
				throw new AssertionError("Wrong NN score trace header.");
			}
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				fields.set(line);
				if(fields.terms()!=5 || fields.parseInt(0)!=++event || fields.parseInt(1)!=32 ||
						fields.parseInt(2)!=Float.floatToRawIntBits(0.8f) || fields.parseFloat(3)!=0.8f || fields.parseInt(4)!=accepted){
					throw new AssertionError("NN trace differs from constant model and exact cutoff decision.");
				}
			}
			if(event!=226){throw new AssertionError("Candidate rows omitted or virtual controls used as real joins: "+event);}
		}finally{if(input.close()){throw new IllegalStateException("Cannot read score trace.");}}
		System.err.println("SAME_K_MODEL_GATE_TRACE_PASS candidates=226 accepted="+accepted);
	}

	/** Writes a small model through native I/O and observes terminal writer errors. */
	private static void write(final String path, final ByteBuilder data){
		final ByteStreamWriter out=new ByteStreamWriter(path, false, false, false);
		out.start(); out.print(data);
		if(out.poisonAndWait()){throw new IllegalStateException("Cannot write model fixture.");}
	}
}
