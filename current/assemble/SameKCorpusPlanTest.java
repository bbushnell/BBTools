package assemble;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;

/** Independent small-manifest checks for fold separation and original K membership. @author Fischl */
public final class SameKCorpusPlanTest {

	/** Checks one planned corpus against its original records without reading sequence. */
	public static void main(final String[] args){
		if(args.length!=2){throw new IllegalArgumentException("Expected source.tsv planned.tsv");}
		final HashMap<String, byte[]> source=new HashMap<String, byte[]>();
		final LineParser1 fields=new LineParser1('\t');
		final ByteFile original=ByteFile.makeByteFile(args[0], true);
		try{
			original.nextLine();
			for(byte[] line=original.nextLine(); line!=null; line=original.nextLine()){
				fields.set(line);
				if(source.put(fields.parseString(0), line)!=null){throw new AssertionError("Duplicate fixture ID");}
			}
		}finally{
			if(original.close()){throw new AssertionError("Source read failure");}
		}
		final HashMap<String, String> groupFold=new HashMap<String, String>();
		final HashSet<String> seen=new HashSet<String>();
		final HashSet<String> folds=new HashSet<String>();
		final ByteFile planned=ByteFile.makeByteFile(args[1], true);
		int rows=0;
		try{
			final byte[] metadata=planned.nextLine(), header=planned.nextLine();
			if(metadata==null || header==null){throw new AssertionError("Missing planned headers");}
			for(byte[] line=planned.nextLine(); line!=null; line=planned.nextLine()){
				fields.set(line);
				if(fields.terms()!=16){throw new AssertionError("Planned row does not have sixteen columns");}
				final String id=fields.parseString(0), group=fields.parseString(3), fold=fields.parseString(12);
				if(!seen.add(id) || !source.containsKey(id)){throw new AssertionError("Lost or duplicated source identity: "+id);}
				int tabs=0, end=0;
				for(; end<line.length; end++){
					if(line[end]=='\t' && ++tabs==12){break;}
				}
				if(!Arrays.equals(source.get(id), Arrays.copyOf(line, end))){throw new AssertionError("Changed original provenance: "+id);}
				final String previous=groupFold.put(group, fold);
				if(previous!=null && !previous.equals(fold)){throw new AssertionError("Group leaked across folds: "+group);}
				if(!fold.equals("train") && !fold.equals("validation") && !fold.equals("test")){
					throw new AssertionError("Unknown fold: "+fold);
				}
				folds.add(fold);
				final int k=fields.parseInt(13);
				boolean found=false;
				for(String token : fields.parseString(9).split(",")){found|=Integer.parseInt(token)==k;}
				if(!found || k<24 || k>127){throw new AssertionError("K outside original schedule: "+id);}
				if(!fields.parseString(14).endsWith("/"+id+"/reads.fq.gz") ||
						!fields.parseString(15).endsWith("/"+id+"/reference.fa")){
					throw new AssertionError("Input paths disagree with genome identity: "+id);
				}
				rows++;
			}
		}finally{
			if(planned.close()){throw new AssertionError("Planned manifest read failure");}
		}
		if(rows!=source.size() || folds.size()!=3){throw new AssertionError("Incomplete rows or empty fold");}
		System.err.println("SAME_K_CORPUS_TEST_PASS rows="+rows+" groups="+groupFold.size());
	}
}
