package stego;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FASTQ;
import stream.Read;
import stream.SamLine;
import stream.SamReadInputStream;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;

/** Encodes arbitrary files in synthetic sequence and recovers them without headers.
 * @author Brian Bushnell, Yaoyao */
public class StegoFast{

	public static void main(String[] args) throws IOException{
		Timer timer=new Timer();
		boolean force=FASTQ.FORCE_INTERLEAVED, test=FASTQ.TEST_INTERLEAVED;
		StegoFast tool=null;
		try{tool=new StegoFast(args); tool.process(timer);}finally{
			FASTQ.FORCE_INTERLEAVED=force;
			FASTQ.TEST_INTERLEAVED=test;
			if(tool!=null){Shared.closeStream(tool.outstream);}
		}
	}

	public StegoFast(String[] args){
		PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : pp.args){
			String[] split=Parse.splitOnFirst(arg, '=');
			String a=split[0].toLowerCase(java.util.Locale.ROOT);
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("in1")){
				if(b==null || b.isEmpty()){throw new IllegalArgumentException("in requires a filename");}
				for(String s : b.split(",", -1)){inputs.add(s);}
			}else if(a.equals("decode")){decode=Parse.parseBoolean(b);
			}else if(a.equals("literal")){literal=b;
			}else if(a.equals("name")){literalName=b;
			}else if(a.equals("len") || a.equals("length")){length=Parse.parseIntKMG(b);
			}else if(a.equals("paired")){paired=b.equalsIgnoreCase("auto") ? null : Parse.parseBoolean(b);
			}else if(a.equals("appendto")){appendTo=b;
			}else if(a.equals("appendto2")){appendTo2=b;
			}else if(a.equals("maxfiles")){maxFiles=Integer.parseInt(b);
			}else if(parser.parse(arg, a, b)){
				// Common BBTools arguments.
			}else{throw new IllegalArgumentException("Unknown parameter "+arg);}
		}
		Parser.processQuality();
		out=parser.out1; out2=parser.out2; in2=parser.in2;
		overwrite=parser.overwrite;
		threadsIn=parser.threadsIn<0 ? 1 : parser.threadsIn;
		threadsOut=parser.threadsOut<0 ? 1 : parser.threadsOut;
		if(parser.append){throw new IllegalArgumentException("Use appendto=reads out=combined, not append=t");}
		if(parser.maxReads>=0){throw new IllegalArgumentException("StegoFast requires complete input; reads= is not supported");}
		if(out==null){throw new IllegalArgumentException("out is required");}
		if(length<128 || length>1000000){throw new IllegalArgumentException("len must be 128..1000000");}
		if(maxFiles<1){throw new IllegalArgumentException("maxfiles must be positive");}
		if(decode){
			if(inputs.size()!=1 || literal!=null || appendTo!=null || appendTo2!=null || out2!=null || overwrite){
				throw new IllegalArgumentException("decode requires one sequence input, optional in2, and out=directory; no overwrite");
			}
			// Header-independent decoding must not require adjacent mate names to match.
			FASTQ.FORCE_INTERLEAVED=false;
			FASTQ.TEST_INTERLEAVED=false;
		}else{
			if(inputs.isEmpty() && literal==null){throw new IllegalArgumentException("Supply in=files or literal=text");}
			if(in2!=null){throw new IllegalArgumentException("Use appendto2 for the second existing read file");}
			if(appendTo2!=null && appendTo==null){throw new IllegalArgumentException("appendto2 requires appendto");}
			if(out2!=null && Boolean.FALSE.equals(paired)){throw new IllegalArgumentException("out2 requires paired metadata");}
			if(!Tools.testOutputFiles(overwrite, false, false, out, out2)){
				throw new IllegalArgumentException("Output exists or cannot be written; use ow=t to replace sequence output");
			}
			for(String input : inputs){checkDistinct(input, out); checkDistinct(input, out2);}
			checkDistinct(appendTo, out); checkDistinct(appendTo, out2);
			checkDistinct(appendTo2, out); checkDistinct(appendTo2, out2);
			checkDistinct(out, out2);
		}
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
	}

	private static void checkDistinct(String a, String b){
		if(a==null || b==null){return;}
		Path x=Paths.get(a).toAbsolutePath().normalize(), y=Paths.get(b).toAbsolutePath().normalize();
		try{
			if(x.equals(y) || (Files.exists(x) && Files.exists(y) && Files.isSameFile(x, y))){
				throw new IllegalArgumentException("Input/output paths must be distinct: "+a+", "+b);
			}
		}catch(IOException e){throw new IllegalArgumentException("Cannot validate paths", e);}
	}

	private void process(Timer timer) throws IOException{
		if(decode){decode();}else{encode();}
		timer.stop();
		outstream.println("Metadata reads: "+metadataReads+"; files: "+files+"; elapsed: "+timer);
	}

	private Streamer reader(FileFormat f1, String second, boolean saveHeader){
		FileFormat f2=FileFormat.testInput(second, FileFormat.FASTQ, null, true, true);
		return StreamerFactory.makeStreamer(f1, f2, true, -1, saveHeader, !f1.samOrBam(), threadsIn);
	}

	private void decode() throws IOException{
		try(StegoArchive.Decoder decoder=new StegoArchive.Decoder(Paths.get(out), maxFiles)){
			try{
				// Read mate files independently; record order and header matching are irrelevant.
				decodeInput(inputs.get(0), decoder);
				if(in2!=null){decodeInput(in2, decoder);}
				files=decoder.finish();
			}catch(IOException | RuntimeException e){
				outstream.println("Extraction failed; temporary evidence retained at "+decoder.workDirectory());
				throw e;
			}
		}
	}

	private void decodeInput(String path, StegoArchive.Decoder decoder) throws IOException{
		FileFormat format=FileFormat.testInput(path, FileFormat.FASTQ, null, true, true);
		Streamer st=reader(format, null, false);
		try{
			st.start();
			if(format.samOrBam()){
				ListNum<SamLine> batch;
				while((batch=st.nextLines())!=null){for(SamLine line : batch.list){accept(line.seq, decoder);}}
			}else{
				ListNum<Read> batch;
				while((batch=st.nextList())!=null){
					for(Read r : batch.list){
						accept(r.bases, decoder);
						if(r.mate!=null){accept(r.mate.bases, decoder);}
					}
				}
			}
		}finally{st.close();}
		if(st.errorState()){throw new IOException("Sequence reader failed: "+path);}
	}

	private void accept(byte[] bases, StegoArchive.Decoder decoder) throws IOException{
		StegoCodec.Chunk chunk=StegoCodec.decode(bases);
		if(chunk!=null){decoder.accept(chunk); metadataReads++;}
	}

	private void encode() throws IOException{
		// Preflight all attachments before opening sequence output.
		java.util.HashSet<String> names=new java.util.HashSet<String>();
		for(String input : inputs){
			Path p=Paths.get(input);
			if(!Files.isRegularFile(p) || !Files.isReadable(p)){throw new IOException("Not a readable file: "+input);}
			String name=p.getFileName().toString();
			StegoArchive.validateName(name);
			if(!names.add(name)){throw new IOException("Duplicate attachment basename: "+name);}
		}
		if(literal!=null){
			StegoArchive.validateName(literalName);
			if(!names.add(literalName)){throw new IOException("Duplicate literal filename: "+literalName);}
		}
		FileFormat f1=FileFormat.testOutput(out, FileFormat.FASTQ, null, true, overwrite, false, true);
		FileFormat f2=FileFormat.testOutput(out2, FileFormat.FASTQ, null, true, overwrite, false, true);
		samOutput=f1.samOrBam();
		if(f2!=null && f2.samOrBam()!=samOutput){throw new IOException("Both mate outputs must use sequence formats or both use alignment formats");}
		FileFormat appendFormat=appendTo==null ? null : FileFormat.testInput(appendTo, FileFormat.FASTQ, null, true, true);
		boolean samInput=appendFormat!=null && appendFormat.samOrBam();
		if(samInput && appendTo2!=null){throw new IOException("Use one SAM/BAM append input; appendto2 is for separate sequence mate files");}
		boolean saveHeader=samInput && f1.samOrBam();
		Streamer st=appendTo==null ? null : reader(appendFormat, appendTo2, saveHeader);
		Writer writer=null;
		boolean complete=false;
		try{
			if(st!=null){st.start();}
			ListNum<SamLine> firstSam=null;
			if(st!=null && samInput){
				firstSam=st.nextLines();// Header publication precedes the first batch or normal EOF.
				if(st.errorState()){throw new IOException("Append input reader failed before first batch");}
			}
			boolean pair=(paired!=null ? paired : out2!=null || (st!=null && st.paired()));
			if(out2!=null && st!=null && !st.paired()){throw new IOException("Separate mate outputs require paired append input");}
			if(st!=null && !samInput && paired!=null && pair!=st.paired()){
				throw new IOException("paired setting conflicts with append input");
			}
			ArrayList<byte[]> headers=saveHeader ? unsortedHeader(SamReadInputStream.getSharedHeader(false)) : null;
			writer=WriterFactory.makeWriter(f1, f2, threadsOut, headers, false);
			writer.start();
			if(st!=null){
				if(samInput){
					for(ListNum<SamLine> batch=firstSam; batch!=null; batch=st.nextLines()){
						if(paired==null && !pair){
							for(SamLine line : batch.list){if(line.hasMate()){pair=true; break;}}
						}
						writer.addLines(new ListNum<SamLine>(batch.list, batchID++));
					}
				}else{
					ListNum<Read> batch;
					while((batch=st.nextList())!=null){writer.add(batch.list, batchID++);}
				}
				st.close();
				if(st.errorState()){throw new IOException("Append input reader failed");}
			}
			for(String input : inputs){
				Path p=Paths.get(input);
				try(InputStream raw=Files.newInputStream(p)){
					encodeFile(p.getFileName().toString(), Files.size(p), raw, pair, writer);
				}
			}
			if(literal!=null){
				byte[] bytes=literal.getBytes(StandardCharsets.UTF_8);
				encodeFile(literalName, bytes.length, new ByteArrayInputStream(bytes), pair, writer);
			}
			flush(writer);
			if(writer.poisonAndWait()){throw new IOException("Sequence writer failed; output may be incomplete");}
			complete=true;
		}finally{
			try{if(!complete && writer!=null){writer.finishError();}}
			finally{if(st!=null){st.close();}}
		}
	}

	/** Preserve dictionaries and other header records without claiming a retained sort. */
	private static ArrayList<byte[]> unsortedHeader(ArrayList<byte[]> original){
		ArrayList<byte[]> copy=new ArrayList<byte[]>();
		boolean hd=false;
		if(original!=null){
			for(byte[] line : original){
				String text=new String(line, StandardCharsets.UTF_8);
				if(text.equals("@HD") || text.startsWith("@HD\t")){
					StringBuilder sb=new StringBuilder("@HD");
					for(String field : text.split("\t")){
						if(!field.equals("@HD") && !field.startsWith("SO:") && !field.startsWith("SS:") && !field.startsWith("GO:")){
							sb.append('\t').append(field);
						}
					}
					copy.add(sb.append("\tSO:unsorted").toString().getBytes(StandardCharsets.UTF_8));
					hd=true;
				}else{copy.add(line.clone());}
			}
		}
		if(!hd){copy.add(0, "@HD\tVN:1.6\tSO:unsorted".getBytes(StandardCharsets.US_ASCII));}
		return copy;
	}

	private void encodeFile(String name, long size, InputStream raw, boolean pair, Writer writer) throws IOException{
		long group;
		do{group=random.nextLong()&Long.MAX_VALUE;}while(usedGroups.contains(Long.valueOf(group)));
		usedGroups.add(Long.valueOf(group));// One identifier allocation per attachment, not per read.
		try(StegoArchive.Payload payload=new StegoArchive.Payload(name, size, raw, length, group)){
			Read pending=null;
			long firstPart=0;
			byte[] buffer=new byte[length/4];
			for(long part=0; part<payload.parts; part++){
				int cap=StegoCodec.capacity(length, group, part), n=0;
				while(n<cap){int read=payload.read(buffer, n, cap-n); if(read<0){break;} n+=read;}
				byte[] bases=StegoCodec.encode(length, group, part, buffer, n);
				Read r=makeRead(bases, null);
				if(!pair){
					r.id=stem(group, part, part);
					add(r, writer);
				}else if(pending==null){pending=r; firstPart=part;
				}else{link(pending, r, stem(group, firstPart, part)); add(pending, writer); pending=null;}
			}
			if(payload.read()!=-1){throw new IOException("Internal metadata part-count mismatch");}
			if(pending!=null){
				Read duplicate=makeRead(pending.bases.clone(), null);
				link(pending, duplicate, stem(group, firstPart, firstPart));
				add(pending, writer);
			}
		}
		files++;
	}

	private static String stem(long group, long first, long second){
		return "metadata_group_"+Long.toHexString(group)+"_part_"+first+(first==second ? "" : "-"+second);
	}
	private Read makeRead(byte[] bases, String name){
		byte[] quality=new byte[bases.length];
		Arrays.fill(quality, (byte)30);
		Read read=new Read(bases, quality, name, numericID++);
		read.setSynthetic(true);
		metadataReads++;
		return read;
	}
	private void link(Read a, Read b, String name){
		a.id=samOutput ? name : name+" 1:";
		b.id=samOutput ? name : name+" 2:";
		a.mate=b; b.mate=a;
		a.setPairnum(0); b.setPairnum(1);
	}
	private void add(Read read, Writer writer){
		batch.add(read);
		if(batch.size()>=200){flush(writer);}
	}
	private void flush(Writer writer){
		if(batch.isEmpty()){return;}
		writer.add(batch, batchID++);
		batch=new ArrayList<Read>(200);
	}

	private final ArrayList<String> inputs=new ArrayList<String>();
	private final ArrayList<Long> usedGroups=new ArrayList<Long>();
	private final SecureRandom random=new SecureRandom();
	private ArrayList<Read> batch=new ArrayList<Read>(200);
	private PrintStream outstream=System.err;
	private String out, out2, in2, appendTo, appendTo2, literal;
	private String literalName="metadata.txt";
	private Boolean paired=null;
	private boolean decode, overwrite, samOutput;
	private int length=150, threadsIn=1, threadsOut=1, maxFiles=128;
	private long batchID, numericID, metadataReads, files;
}
