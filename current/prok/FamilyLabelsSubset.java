package prok;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.PreParser;
import prot.ProteinFamilyAssigner;
import prot.ProteinSequence;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Transfers existing labels only to an exact ordered subset of unchanged proteins.
 * Validates both model packages and every selected consensus before rebinding
 * local ranks and decoded roster/consensus pins. All other label fields survive.
 * This does not infer names for newly trained or changed sequences.
 * @author Keqing
 */
public final class FamilyLabelsSubset {

	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, FamilyLabelsSubset.class, false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String[] args) throws IOException{
		Path parent=null, subset=null, labels=null, output=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			require(eq>0 && eq<arg.length()-1, "Expected key=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT);
			final Path value=Paths.get(arg.substring(eq+1));
			if(key.equals("parent")){require(parent==null, "Duplicate parent"); parent=value;}
			else if(key.equals("models")){require(subset==null, "Duplicate models"); subset=value;}
			else if(key.equals("labels")){require(labels==null, "Duplicate labels"); labels=value;}
			else if(key.equals("out")){require(output==null, "Duplicate out"); output=value;}
			else{throw new IllegalArgumentException("Unknown label subset parameter: "+key);}
		}
		require(parent!=null && subset!=null && labels!=null && output!=null,
			"Required: parent=MODELS models=SUBSET_MODELS labels=PARENT_LABELS out=NEW_RAW_TSV");
		validateOutput(output);
		final Path target=output.toFile().getCanonicalFile().toPath();
		require(!target.startsWith(parent.toRealPath()) && !target.startsWith(subset.toRealPath()),
			"Label output must not modify an immutable model directory");
		Shared.setThreads(1); Shared.AMINO_IN=true;
		final ProteinFamilyAssigner original=ProteinFamilyAssigner.load(parent);
		final ProteinFamilyAssigner selected=ProteinFamilyAssigner.load(subset);
		FamilyLabels.load(labels.toString(), original);
		final byte[][] a=consensus(parent, original), b=consensus(subset, selected);
		final int[] local=new int[original.familyCount()]; Arrays.fill(local, -1);
		int next=0;
		for(int r=0; r<selected.familyCount(); r++){
			while(next<original.familyCount() && original.familyId(next)!=selected.familyId(r)){next++;}
			require(next<original.familyCount(), "Subset family is absent or out of parent order: "+selected.familyId(r));
			require(original.representative(next).equals(selected.representative(r))
				&& original.memberCount(next)==selected.memberCount(r) && Arrays.equals(a[next], b[r]),
				"Labels cannot transfer to changed family/consensus: "+selected.familyId(r));
			local[next++]=r;
		}
		final int written=write(labels, output, selected.familyCount(), selected.identityPin(), selected.consensusPin(), local);
		FamilyLabels.load(output.toString(), selected);
		System.err.println("FAMILY_LABELS_SUBSET_PASS families="+selected.familyCount()+" labels="+written
			+" unchanged_sequences=true retained_label_fields=true");
	}

	/** Output must be reopenable for validation; special streams and codecs are unsupported. */
	public static void validateOutput(final Path output){
		final String name=output.toString();
		// ReadWrite.getRawOutputStream redirects stderr and rewrites '|', beyond FileFormat's classification.
		require(!name.equals("stderr") && !name.startsWith("stderr.") && name.indexOf('|')<0,
			"Output must be an ordinary reopenable filename: "+output);
		final FileFormat format=FileFormat.testOutput(name, FileFormat.TEXT, null, false, false, false, false);
		require(format.file() && !format.compressed() && !Files.exists(output),
			"Output must be a fresh ordinary uncompressed file; compress it when packaging: "+output);
	}

	/** Uses the existing FASTA streamer; the library provides the expected roster. */
	private static byte[][] consensus(final Path directory, final ProteinFamilyAssigner library){
		final Path compressed=directory.resolve("consensus.faa.gz");
		final Path path=Files.isRegularFile(compressed) ? compressed : directory.resolve("consensus.faa");
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(path.toString(), FileFormat.FASTA,
			null, true, false), null, true, -1);
		final byte[][] sequences=new byte[library.familyCount()][]; int rank=0;
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null && batch.size()>0; batch=input.nextList()){
				for(Read read : batch){
					require(read.mate==null && read.id!=null && read.bases!=null, "Expected named unpaired protein FASTA records");
					int end=0; while(end<read.id.length() && !Character.isWhitespace(read.id.charAt(end))){end++;}
					require(rank<sequences.length && library.representative(rank).equals(read.id.substring(0, end)),
						"Consensus differs from validated model order at rank "+rank);
					sequences[rank]=new ProteinSequence(library.representative(rank), read.bases).enc; rank++;
				}
			}
		}finally{input.close();}
		require(!input.errorState() && rank==sequences.length, "Consensus read failed or family count differs");
		return sequences;
	}

	/** Rewrites only two binding headers and active_index, preserving all other bytes. */
	public static int write(final Path labels, final Path output, final int families, final String identityPin,
			final String consensusPin, final int[] local){
		require(families>0 && local!=null, "Label subset writer needs a nonempty validated mapping");
		final ByteFile reader=ByteFile.makeByteFile(labels.toString(), false);
		final java.util.ArrayList<byte[]> headers=new java.util.ArrayList<byte[]>();
		final byte[][] rows=new byte[families][];
		final LineParser1 lp=new LineParser1('\t'); final ByteBuilder row=new ByteBuilder();
		int rankColumn=-1, width=0, written=0;
		try{
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				if(line.length==0){continue;} lp.set(line);
				if(line[0]=='#'){
					if(lp.termEquals("#family_identity_sha80", 0)){
						headers.add(new ByteBuilder().append("#family_identity_sha80\t").append(identityPin).toBytes());
					}else if(lp.termEquals("#consensus_sha80", 0)){
						headers.add(new ByteBuilder().append("#consensus_sha80\t").append(consensusPin).toBytes());
					}else{headers.add(line);}
					continue;
				}
				if(width==0){
					width=lp.terms(); headers.add(line);
					for(int c=0; c<width; c++){if(lp.termEquals("active_index", c)){rankColumn=c;}}
					require(rankColumn>=0, "Validated labels lost their active_index column"); continue;
				}
				require(lp.terms()==width, "Label width changed after validation");
				final int source=lp.parseInt(rankColumn);
				require(source>=0 && source<local.length, "Label rank changed after validation");
				final int destination=local[source]; if(destination<0){continue;}
				require(rows[destination]==null, "Duplicate subset label rank: "+destination);
				row.clear();
				for(int c=0; c<width; c++){
					if(c>0){row.tab();}
					if(c==rankColumn){row.append(destination);}else{lp.appendTerm(row, c);}
				}
				rows[destination]=row.toBytes(); written++;
			}
		}finally{require(!reader.close(), "Label input read failed");}
		require(width>0, "Missing label header");
		final ByteStreamWriter writer=ByteStreamWriter.makeBSW(output.toString(), false, false, true);
		try{
			for(byte[] header : headers){writer.println(header);}
			for(byte[] record : rows){if(record!=null){writer.println(record);}}
		}finally{require(!writer.poisonAndWait(), "Subset label write failed");}
		return written;
	}

	private static void require(final boolean ok, final String reason){if(!ok){throw new IllegalArgumentException(reason);}}
}
