package test.stego;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

import stego.StegoArchive;
import stego.StegoCodec;
import fileIO.FileFormat;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Round-trip and adversarial regression tests; runs in an isolated directory.
 * @author Yaoyao */
public class StegoFastTest{

	public static void main(String[] args) throws Exception{
		StegoFastTest test=new StegoFastTest(Paths.get(args[0]).toAbsolutePath());
		System.out.println("Artifacts: "+test.work);
		test.codecTests();
		test.archiveTests();
		System.out.println("Codec and archive checks passed: "+test.checks);
		test.integrationTests();
		System.out.println("PASS: "+test.checks+" checks; artifacts: "+test.work);
	}
	private StegoFastTest(Path root) throws IOException{
		repo=root;
		work=Files.createTempDirectory(repo.getParent(), "stegofast-test-");
	}
	private void check(boolean condition, String message){
		checks++;
		if(!condition){throw new AssertionError(message);}
	}
	private interface Action{void run() throws Exception;}
	private void rejects(Action action, String message) throws Exception{
		boolean rejected=false;
		try{action.run();}catch(IOException | IllegalArgumentException | ArithmeticException e){rejected=true;}
		check(rejected, message);
	}
	private static byte[] reverse(byte[] source){
		byte[] b=new byte[source.length];
		for(int i=0; i<b.length; i++){
			switch(source[source.length-1-i]){
				case 'A': b[i]='T'; break;
				case 'C': b[i]='G'; break;
				case 'G': b[i]='C'; break;
				case 'T': b[i]='A'; break;
				default: throw new AssertionError("Unexpected base");
			}
		}
		return b;
	}
	private void codecTests() throws Exception{
		check(StegoCodec.decode("ACGTACGTACGT".getBytes(StandardCharsets.US_ASCII))==null, "Ordinary read");
		for(long value : new long[]{0,1,3,4,15,16,255,256,65535,65536,Long.MAX_VALUE}){
			byte[] payload={0,27,(byte)255};
			byte[] bases=StegoCodec.encode(151, value, value, payload, payload.length);
			StegoCodec.Chunk c=StegoCodec.decode(bases);
			check(c.group==value && c.part==value && Arrays.equals(payload, c.data), "Integer boundary "+value);
			check(Arrays.equals(payload, StegoCodec.decode(reverse(bases)).data), "Reverse complement");
			check(Arrays.equals(payload, StegoCodec.decode(new String(bases, StandardCharsets.US_ASCII)
				.toLowerCase(java.util.Locale.ROOT).getBytes(StandardCharsets.US_ASCII)).data), "Lowercase");
		}
		byte[] known=StegoCodec.encode(150, 1, 0, new byte[]{27}, 1);
		String prefix=StegoCodec.MAGIC+"AC"+"AAAC"+"AAAA"+"AAAC"+"ACGT";
		check(new String(known, StandardCharsets.US_ASCII).startsWith(prefix), "Hand-derived framing/base mapping");
		for(int i=20; i<known.length; i++){
			final byte[] damaged=known.clone();
			damaged[i]=(byte)(damaged[i]=='A' ? 'C' : 'A');
			rejects(()->StegoCodec.decode(damaged), "Every nonmagic position is protected: "+i);
		}
		rejects(()->StegoCodec.decode(Arrays.copyOf(known, known.length-1)), "Truncation");
	}
	private ArrayList<byte[]> encoded(byte[] bytes, String name, int length, long group) throws IOException{
		ArrayList<byte[]> reads=new ArrayList<byte[]>();
		try(StegoArchive.Payload in=new StegoArchive.Payload(name, bytes.length, new ByteArrayInputStream(bytes), length, group)){
			for(long part=0; part<in.parts; part++){
				byte[] buffer=new byte[StegoCodec.capacity(length, group, part)];
				int used=0;
				while(used<buffer.length){int n=in.read(buffer, used, buffer.length-used); if(n<0){break;} used+=n;}
				reads.add(StegoCodec.encode(length, group, part, buffer, used));
			}
			check(in.read()==-1, "Part count fits manifest stream");
		}
		return reads;
	}
	private void extract(List<byte[]> reads, Path dir) throws IOException{
		try(StegoArchive.Decoder decoder=new StegoArchive.Decoder(dir, 128)){
			for(byte[] b : reads){decoder.accept(StegoCodec.decode(b));}
			check(decoder.finish()==1, "One recovered file");
		}
	}
	private void archiveTests() throws Exception{
		Random rng=new Random(67931);
		for(int length : new int[]{128,150,151,256,400}){
			for(int size : new int[]{0,1,7,16,255,4096}){
				byte[] data=new byte[size]; rng.nextBytes(data);
				String name="notes-\u96e8.bin";
				ArrayList<byte[]> reads=encoded(data, name, length, Long.MAX_VALUE);
				reads.add(reads.get(0).clone());
				Collections.shuffle(reads, rng);
				for(int i=0; i<reads.size(); i+=2){reads.set(i, reverse(reads.get(i)));}
				Path dir=work.resolve("unit-"+length+"-"+size);
				extract(reads, dir);
				check(Arrays.equals(data, Files.readAllBytes(dir.resolve(name))), "Shuffled bytes "+size+" at "+length);
			}
		}
		byte[] data=new byte[300]; rng.nextBytes(data);
		ArrayList<byte[]> reads=encoded(data, "original.bin", 150, 2718);
		ArrayList<byte[]> missing=new ArrayList<byte[]>(reads); missing.remove(2);
		rejects(()->extract(missing, work.resolve("missing")), "Missing middle part");
		ArrayList<byte[]> shortFile=new ArrayList<byte[]>(reads); shortFile.remove(shortFile.size()-1);
		rejects(()->extract(shortFile, work.resolve("missing-last")), "Missing last part");
		StegoCodec.Chunk original=StegoCodec.decode(reads.get(3));
		original.data[0]^=1;
		ArrayList<byte[]> conflicting=new ArrayList<byte[]>(reads);
		conflicting.add(StegoCodec.encode(150, 2718, 3, original.data, original.data.length));
		rejects(()->extract(conflicting, work.resolve("conflict")), "Conflicting valid duplicate");
		ArrayList<byte[]> changed=new ArrayList<byte[]>(reads);
		StegoCodec.Chunk last=StegoCodec.decode(changed.get(changed.size()-1));
		last.data[last.data.length-1]^=1;
		changed.set(changed.size()-1, StegoCodec.encode(150, 2718, last.part, last.data, last.data.length));
		rejects(()->extract(changed, work.resolve("file-crc")), "File CRC despite valid read CRCs");
		Path existing=work.resolve("existing"); Files.createDirectory(existing);
		Files.write(existing.resolve("original.bin"), new byte[]{7});
		rejects(()->extract(reads, existing), "Existing destination");
		check(Arrays.equals(new byte[]{7}, Files.readAllBytes(existing.resolve("original.bin"))), "Existing file untouched");
		for(String name : new String[]{"../escape","/absolute","..","x\\y","C:foo"}){
			rejects(()->StegoArchive.validateName(name), "Unsafe basename");
		}
	}
	private void command(boolean success, String script, String... args) throws Exception{
		ArrayList<String> command=new ArrayList<String>();
		command.add("bash"); command.add(repo.resolve(script).toString());
		command.add("-Xmx256m"); command.add("t=2");
		command.addAll(Arrays.asList(args));
		Path log=work.resolve("command-"+(commands++)+".log");
		Process p=new ProcessBuilder(command).directory(repo.toFile()).redirectErrorStream(true).redirectOutput(log.toFile()).start();
		if(!p.waitFor(90, TimeUnit.SECONDS)){p.destroyForcibly(); p.waitFor(); throw new AssertionError("Hung: "+command+"; "+log);}
		check((p.exitValue()==0)==success, "Exit "+p.exitValue()+": "+command+"; log="+log);
	}
	private void integrationTests() throws Exception{
		byte[] binary=new byte[4096]; new Random(42).nextBytes(binary);
		Path original=work.resolve("binary.bin"); Files.write(original, binary);
		Path empty=work.resolve("empty.bin"); Files.write(empty, new byte[0]);
		ByteArrayOutputStream compressed=new ByteArrayOutputStream();
		try(GZIPOutputStream gzip=new GZIPOutputStream(compressed)){gzip.write(binary);}
		Path gz=work.resolve("opaque.gz"); Files.write(gz, compressed.toByteArray());
		for(String ext : new String[]{"fq.gz","fa","sam","bam"}){
			Path encoded=work.resolve("encoded."+ext), recovered=work.resolve("decoded-"+ext);
			command(true, "stegofast.sh", "in="+original+","+empty+","+gz, "literal=a=b \u96e8", "name=message.txt",
				"out="+encoded, "paired=t");
			command(true, "stegofast.sh", "in="+encoded, "decode", "out="+recovered);
			if(ext.equals("sam") || ext.equals("bam")){checkAlignmentNames(encoded);}
			check(Arrays.equals(binary, Files.readAllBytes(recovered.resolve("binary.bin"))), "Binary through "+ext);
			check(Files.size(recovered.resolve("empty.bin"))==0, "Empty through "+ext);
			check(Arrays.equals(compressed.toByteArray(), Files.readAllBytes(recovered.resolve("opaque.gz"))), "Opaque gzip through "+ext);
			check(new String(Files.readAllBytes(recovered.resolve("message.txt")), StandardCharsets.UTF_8).equals("a=b \u96e8"), "Literal through "+ext);
		}
		Path paired=work.resolve("paired.fq");
		command(true, "stegofast.sh", "in="+original, "out="+paired, "paired=t");
		command(true, "reformat.sh", "in="+paired, "out="+work.resolve("verified.fq"), "verifyinterleaved=t");
		List<String> fq=Files.readAllLines(paired, StandardCharsets.US_ASCII);
		ArrayList<String> scrambled=new ArrayList<String>();
		for(int i=fq.size()-4; i>=0; i-=4){
			scrambled.add(">SRA_changed_"+i);
			scrambled.add(new String(reverse(fq.get(i+1).getBytes(StandardCharsets.US_ASCII)), StandardCharsets.US_ASCII));
		}
		Path fasta=work.resolve("scrambled.fa"); Files.write(fasta, scrambled, StandardCharsets.US_ASCII);
		command(true, "stegofast.sh", "in="+fasta, "decode", "out="+work.resolve("scrambled-decoded"));
		check(Arrays.equals(binary, Files.readAllBytes(work.resolve("scrambled-decoded/binary.bin"))), "Reorder, headers, qualities, orientation");
		ArrayList<String> mapped=new ArrayList<String>();
		mapped.add("@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural");
		mapped.add("@SQ\tSN:synthetic\tLN:100000");
		for(int i=fq.size()-4; i>=0; i-=4){
			String sequence=new String(reverse(fq.get(i+1).getBytes(StandardCharsets.US_ASCII)), StandardCharsets.US_ASCII);
			mapped.add("rewritten_"+i+"\t16\tsynthetic\t"+(fq.size()-i)+"\t30\t5S140M5S\t*\t0\t0\t"+sequence+"\t"+fq.get(i+3));
		}
		Path mappedSam=work.resolve("mapped.sam"); Files.write(mappedSam, mapped, StandardCharsets.US_ASCII);
		command(true, "stegofast.sh", "in="+mappedSam, "decode", "out="+work.resolve("mapped-sam-decoded"));
		check(Arrays.equals(binary, Files.readAllBytes(work.resolve("mapped-sam-decoded/binary.bin"))), "Mapped SAM without MD tag");
		command(true, "stegofast.sh", "literal=transport", "appendto="+mappedSam, "out="+work.resolve("mapped.bam"));
		command(true, "stegofast.sh", "in="+work.resolve("mapped.bam"), "decode", "out="+work.resolve("mapped-decoded"));
		check(Arrays.equals(binary, Files.readAllBytes(work.resolve("mapped-decoded/binary.bin"))), "Mapped reverse-strand soft-clipped BAM");
		Path copiedSam=work.resolve("appended-mapped.sam");
		command(true, "stegofast.sh", "literal=extra", "appendto="+mappedSam, "out="+copiedSam);
		List<String> copied=Files.readAllLines(copiedSam, StandardCharsets.US_ASCII);
		check(copied.get(0).equals("@HD\tVN:1.6\tSO:unsorted"), "Appended sort declaration");
		check(copied.get(1).equals(mapped.get(1)), "Reference dictionary retained");
		check(copied.subList(2, mapped.size()).equals(mapped.subList(2, mapped.size())), "Mapped record fields unchanged");
		Path r1=work.resolve("r1.fq"), r2=work.resolve("r2.fq");
		command(true, "stegofast.sh", "literal=split mates", "out="+r1, "out2="+r2);
		command(true, "stegofast.sh", "in="+r1, "in2="+r2, "decode", "out="+work.resolve("split-decoded"));
		Path reads=work.resolve("ordinary.fq");
		Files.write(reads, Arrays.asList("@biological", "ACGTACGT", "+", "????????"), StandardCharsets.US_ASCII);
		Path combined=work.resolve("combined.fq");
		command(true, "stegofast.sh", "literal=embedded", "appendto="+reads, "out="+combined);
		List<String> combinedLines=Files.readAllLines(combined, StandardCharsets.US_ASCII);
		check(combinedLines.subList(0,4).equals(Files.readAllLines(reads, StandardCharsets.US_ASCII)), "Original read unchanged");
		command(true, "stegofast.sh", "in="+combined, "decode", "out="+work.resolve("append-decoded"));
		command(false, "stegofast.sh", "literal=oops", "appendto="+reads, "out="+reads, "ow=t");
		command(false, "stegofast.sh", "in="+combined, "decode", "out="+work.resolve("append-decoded"));
		command(true, "reformat.sh", "in="+paired, "out="+work.resolve("mapped-container.bam"));
		command(true, "stegofast.sh", "literal=more", "name=extra.txt", "appendto="+work.resolve("mapped-container.bam"), "out="+work.resolve("combined.bam"));
		command(true, "stegofast.sh", "in="+work.resolve("combined.bam"), "decode", "out="+work.resolve("bam-append-decoded"));
		check(Arrays.equals(binary, Files.readAllBytes(work.resolve("bam-append-decoded/binary.bin"))), "BAM append preserves earlier group");
	}
	private void checkAlignmentNames(Path input) throws IOException{
		FileFormat format=FileFormat.testInput(input.toString(), FileFormat.SAM, null, true, true);
		Streamer st=StreamerFactory.makeStreamer(format, null, true, -1, false, false, 1);
		String previous=null;
		int count=0;
		try{
			st.start();
			ListNum<SamLine> batch;
			while((batch=st.nextLines())!=null){
				for(SamLine line : batch.list){
					check(line.qname.indexOf(' ')<0 && line.qname.indexOf('\t')<0, "Whitespace-free alignment QNAME");
					check(line.hasMate() && !line.mapped(), "Synthetic unmapped paired flags");
					if((count&1)==0){previous=line.qname; check((line.flag&64)!=0, "First-mate flag");}
					else{check(previous.equals(line.qname) && (line.flag&128)!=0, "Matching mate QNAME and second flag");}
					count++;
				}
			}
		}finally{st.close();}
		check(!st.errorState() && count>0 && (count&1)==0, "Complete alignment pairs");
	}
	private final Path repo, work;
	private int checks, commands;
}
