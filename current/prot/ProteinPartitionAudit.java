package prot;

import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Verifies round-robin FASTA partitions against every original header and residue.
 * Reconstructs source order through native sequence streams without retaining the
 * corpus. Amino validation is disabled so malformed source symbols cannot be
 * silently normalized by the verification reader. This tool repairs nothing.
 * @author Keqing
 */
public final class ProteinPartitionAudit {
	public static void main(String[] args){
		try{
			Shared.setThreads(1); Shared.AMINO_IN=true; Shared.TRIM_READ_DESCRIPTION=false; Read.VALIDATE_IN_CONSTRUCTOR=false;
			PreParser pp=new PreParser(args, ProteinPartitionAudit.class, false);
			try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable error){error.printStackTrace(); System.exit(1);}
	}
	private static void run(HashMap<String,String> o){
		for(String key : o.keySet()){require(Arrays.asList("in", "parts", "ways", "expected", "out").contains(key), "Unknown partition audit argument: "+key);}
		final int ways=Integer.parseInt(HmmComparisonData.required(o, "ways"));
		final long expected=Long.parseLong(HmmComparisonData.required(o, "expected"));
		final String pattern=HmmComparisonData.required(o, "parts");
		require(ways>0 && ways<=64 && expected>0 && pattern.indexOf('%')>=0 && pattern.indexOf('%')==pattern.lastIndexOf('%'),
			"Partition audit needs1..64 ways, positive expected count and exactly one filename placeholder");
		final Cursor source=new Cursor(HmmComparisonData.required(o, "in")); final Cursor[] parts=new Cursor[ways];
		final long[] counts=new long[ways], residues=new long[ways]; long total=0;
		try{
			for(int i=0; i<ways; i++){parts[i]=new Cursor(pattern.replace("%", Integer.toString(i)));}
			int part=0;
			for(Read original=source.next(); original!=null; original=source.next()){
				final Read copy=parts[part].next();
				if(copy==null || !original.id.equals(copy.id) || !Arrays.equals(original.bases, copy.bases)){
					throw new IllegalArgumentException("Partition header/residue mismatch at source record "+total+" in partition "+part);
				}
				require(original.mate==null && copy.mate==null && original.bases!=null && original.bases.length>0, "Partition audit expects nonempty unpaired proteins");
				counts[part]++; residues[part]+=original.bases.length; total++;
				part++; if(part==ways){part=0;}
				if(total%1000000==0){System.err.println("PROTEIN_PARTITION_AUDIT_PROGRESS records="+total);}
			}
			require(total==expected, "Source count differs from the partition contract: expected="+expected+" actual="+total);
			for(int i=0; i<ways; i++){require(parts[i].next()==null, "Extra records remain in partition "+i);}
		}finally{
			source.close(); for(Cursor part : parts){if(part!=null){part.close();}}
		}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(o, "out"));
		out.println("partition\trecords\traw_residues"); final ByteBuilder row=new ByteBuilder();
		for(int i=0; i<ways; i++){out.println(row.clear().append(i).tab().append(counts[i]).tab().append(residues[i]));}
		HmmComparisonData.close(out);
		System.err.println("PROTEIN_PARTITION_AUDIT_PASS records="+total+" partitions="+ways+" headers=true residues=true");
	}

	/** Owns one bounded native reader and its current batch. Returned Reads remain caller-owned. */
	private static final class Cursor{
		Cursor(String file){stream=StreamerFactory.makeStreamer(FileFormat.testInput(file, FileFormat.FASTA, null, true, false), null, true, -1); stream.start();}
		Read next(){
			while(list==null || index>=list.size()){
				if(ended){return null;}
				list=stream.nextList(); index=0;
				if(list==null || list.size()==0){ended=true; return null;}
			}
			return list.list.get(index++);
		}
		void close(){stream.close(); require(!stream.errorState(), "Partition input stream failed");}
		final Streamer stream;
		ListNum<Read> list;
		int index;
		boolean ended;
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
}
