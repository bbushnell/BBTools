package bin;

import java.io.PrintStream;
import java.util.ArrayList;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.SamLine;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import structures.ListNum;

/**
 * Renames CAMI benchmark contigs to include tid_ labels.
 * Reads a CAMI binning_gs.tsv key file and appends _tid_TAXID
 * to contig names in both FASTA and SAM/BAM files.
 *
 * Usage: renamecami.sh in=contigs.fa sam=mapped.bam key=binning_gs.tsv
 *        out=renamed.fa outsam=renamed.bam
 *
 * @author Eru, Brian Bushnell
 * @date July 24, 2026
 */
public class RenameCAMI {

	public static void main(String[] args) {
		Timer t=new Timer();
		RenameCAMI x=new RenameCAMI(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public RenameCAMI(String[] args) {

		{
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		Parser parser=new Parser();
		for(int i=0; i<args.length; i++) {
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")) {b=null;}

			if(a.equals("key") || a.equals("cami")) {
				keyFile=b;
			} else if(a.equals("sam") || a.equals("insam") || a.equals("bam")) {
				samIn=b;
			} else if(a.equals("outsam") || a.equals("outbam")) {
				samOut=b;
			} else if(a.equals("dropunlabeled") || a.equals("drop")) {
				dropUnlabeled=Parse.parseBoolean(b);
			} else if(a.equals("verbose")) {
				verbose=Parse.parseBoolean(b);
			} else if(parser.parse(arg, a, b)) {
				//do nothing
			} else {
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		in1=parser.in1;
		out1=parser.out1;
		overwrite=parser.overwrite;

		if(keyFile==null) {
			throw new RuntimeException("Error - a CAMI key file is required (key=binning_gs.tsv).");
		}
	}

	void process(Timer t) {
		keyMap=loadKeyFile(keyFile);
		outstream.println("Loaded "+keyMap.size()+" contig-to-taxid mappings from "+keyFile);

		long fastaRenamed=0, fastaDropped=0;
		long samRenamed=0, samDropped=0;

		if(in1!=null && out1!=null) {
			long[] counts=renameFasta();
			fastaRenamed=counts[0];
			fastaDropped=counts[1];
		}

		if(samIn!=null && samOut!=null) {
			long[] counts=renameSam();
			samRenamed=counts[0];
			samDropped=counts[1];
		}

		t.stop();
		if(in1!=null) {
			outstream.println("FASTA: "+fastaRenamed+" contigs renamed, "+fastaDropped+" dropped.");
		}
		if(samIn!=null) {
			outstream.println("SAM: "+samRenamed+" alignments renamed, "+samDropped+" dropped.");
		}
		outstream.println("Time: "+t);
	}

	/**
	 * Loads a CAMI binning_gs.tsv file.
	 * Format: lines starting with @ are headers/metadata.
	 * Data lines are tab-separated: SEQUENCEID BINID TAXID [_LENGTH]
	 */
	private HashMap<String, Integer> loadKeyFile(String fname) {
		HashMap<String, Integer> map=new HashMap<String, Integer>();
		ByteFile bf=ByteFile.makeByteFile(fname, false);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(line.length<1 || line[0]=='@' || line[0]=='#') {continue;}
			String s=new String(line);
			String[] tabs=s.split("\t");
			if(tabs.length<3) {continue;}
			String contigName=tabs[0].trim();
			int taxid=-1;
			try {
				taxid=Integer.parseInt(tabs[2].trim());
			} catch(NumberFormatException e) {
				continue;
			}
			if(taxid>0) {map.put(contigName, taxid);}
		}
		bf.close();
		return map;
	}

	private long[] renameFasta() {
		FileFormat ffin=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
		FileFormat ffout=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, false, true);

		Streamer st=StreamerFactory.makeStreamer(ffin, null, true, -1, false, true, -1);
		Writer fw=WriterFactory.makeWriter(ffout, null, -1, null, false);

		st.start();
		fw.start();

		long renamed=0, dropped=0;
		for(ListNum<Read> ln=st.nextList(); ln!=null; ln=st.nextList()) {
			ArrayList<Read> list=ln.list;
			for(int i=0; i<list.size(); i++) {
				Read r=list.get(i);
				String name=r.id;
				String baseName=name.split("\\s+")[0];
				Integer tid=keyMap.get(baseName);
				if(tid!=null) {
					r.id=baseName+"_tid_"+tid;
					renamed++;
				} else if(dropUnlabeled) {
					list.set(i, null);
					dropped++;
				}
			}
			if(dropUnlabeled) {
				list.removeIf(r -> r==null);
			}
			fw.addReads(ln);
		}

		fw.poisonAndWait();
		st.close();
		return new long[] {renamed, dropped};
	}

	private long[] renameSam() {
		FileFormat ffin=FileFormat.testInput(samIn, FileFormat.SAM, null, true, true);
		FileFormat ffout=FileFormat.testOutput(samOut, FileFormat.SAM, null, true, overwrite, false, true);

		SamLine.SET_FROM_OK=true;

		Streamer st=StreamerFactory.makeStreamer(ffin, null, true, -1, true, false, -1);
		Writer fw=WriterFactory.makeWriter(ffout, null, -1, null, true);

		st.start();
		fw.start();

		long renamed=0, dropped=0;
		for(ListNum<SamLine> ln=st.nextLines(); ln!=null; ln=st.nextLines()) {
			ArrayList<SamLine> list=ln.list;
			for(int i=0; i<list.size(); i++) {
				SamLine sl=list.get(i);
				String rn=sl.rnameS();
				if(rn!=null && !rn.equals("*")) {
					Integer tid=keyMap.get(rn);
					if(tid!=null) {
						sl.setRname((rn+"_tid_"+tid).getBytes(StandardCharsets.US_ASCII));
						renamed++;
					} else if(dropUnlabeled) {
						list.set(i, null);
						dropped++;
					}
				}
			}
			if(dropUnlabeled) {
				list.removeIf(sl -> sl==null);
			}
			fw.addLines(ln);
		}

		fw.poisonAndWait();
		st.close();
		return new long[] {renamed, dropped};
	}

	private String in1=null;
	private String out1=null;
	private String samIn=null;
	private String samOut=null;
	private String keyFile=null;
	private boolean dropUnlabeled=false;
	private boolean overwrite=false;

	private HashMap<String, Integer> keyMap;
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
}
