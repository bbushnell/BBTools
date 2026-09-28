package prok;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.FileFormat;
import map.LongHashSet;
import map.LongObjectMap;
import parse.LineParser1;

/** Immutable per-seed endpoint-offset resource used by conserved-RNA voting. */
public final class SeedOffsetTable {

	public static SeedOffsetTable load(String path, int expectedK, LongHashSet expectedSeeds){
		return load(path, expectedK, expectedSeeds, null);
	}

	public static SeedOffsetTable load(String path, int expectedK, LongHashSet expectedSeeds,
			String expectedSeedSetSha80){
		if(path==null){throw new IllegalArgumentException("Seed-offset table path is null");}
		if(expectedK<1 || expectedK>31){throw new IllegalArgumentException("Invalid expected seed K: "+expectedK);}
		if(expectedSeeds==null || expectedSeeds.size()<1){throw new IllegalArgumentException("Expected seed set is empty");}
		if(expectedSeedSetSha80!=null && !validSha80(expectedSeedSetSha80)){
			throw new IllegalArgumentException("Invalid expected seed-set sha80: "+expectedSeedSetSha80);
		}

		final ByteFile bf=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		final LineParser1 lp=new LineParser1((byte)'\t');
		final LongObjectMap<KmerInfo> map=new LongObjectMap<KmerInfo>(Math.max(16, expectedSeeds.size()*2), KmerInfo.class);
		int version=-1, k=-1, declaredSeeds=-1;
		long trainingGenes=-1;
		String seedSetSha80=null;
		boolean sawVersion=false, sawK=false, sawSeedCount=false, sawTrainingGenes=false, sawSeedSha=false;
		boolean sawHeader=false;
		byte[] line;
		while((line=bf.nextLine())!=null){
			if(line.length<1){continue;}
			if(Arrays.equals(line, HEADER_BYTES)){
				if(sawHeader){fail("Duplicate data header");}
				sawHeader=true;
				continue;
			}
			lp.set(line);
			if(line[0]=='#'){
				if(sawHeader){fail("Metadata appears after data header: "+ascii(line));}
				if(lp.terms()!=2){fail("Malformed metadata line: "+ascii(line));}
				if(lp.termEquals("#SeedOffsetTable", 0)){
					if(sawVersion){fail("Duplicate metadata key: #SeedOffsetTable");} sawVersion=true; version=lp.parseInt(1);
				}else if(lp.termEquals("#K", 0)){
					if(sawK){fail("Duplicate metadata key: #K");} sawK=true; k=lp.parseInt(1);
				}else if(lp.termEquals("#SeedCount", 0)){
					if(sawSeedCount){fail("Duplicate metadata key: #SeedCount");} sawSeedCount=true; declaredSeeds=lp.parseInt(1);
				}else if(lp.termEquals("#TrainingGenes", 0)){
					if(sawTrainingGenes){fail("Duplicate metadata key: #TrainingGenes");} sawTrainingGenes=true; trainingGenes=lp.parseLong(1);
				}else if(lp.termEquals("#SeedSetSha80", 0)){
					if(sawSeedSha){fail("Duplicate metadata key: #SeedSetSha80");} sawSeedSha=true; seedSetSha80=ascii(lp.parseByteArray(1));
				}
				else{fail("Unknown seed-offset metadata key: "+ascii(lp.parseByteArray(0)));}
				continue;
			}
			if(!sawHeader){fail("Data row appears before exact header");}
			if(lp.terms()!=6){fail("Expected six A48 columns, found "+lp.terms()+": "+ascii(line));}
			final long seed=lp.parseLongA48(0);
			final long left10=lp.parseLongA48(1), right10=lp.parseLongA48(2);
			final long leftSd10=lp.parseLongA48(3), rightSd10=lp.parseLongA48(4), count=lp.parseLongA48(5);
			if(seed<0){fail("Negative packed seed: "+seed);}
			if(left10<0 || right10<0 || leftSd10<0 || rightSd10<0 || count<0){
				fail("Negative offset/dispersion/count for seed "+seed);
			}
			if(count==0 && (left10!=0 || right10!=0 || leftSd10!=0 || rightSd10!=0)){
				fail("Untrained seed "+seed+" has nonzero coordinates");
			}
			if(!expectedSeeds.contains(seed)){fail("Extra seed row absent from expected set: "+seed);}
			final KmerInfo info=new KmerInfo(left10/10.0, right10/10.0, leftSd10/10.0, rightSd10/10.0, count);
			if(map.put(seed, info)!=null){fail("Duplicate seed row: "+seed);}
		}
		if(bf.close()){throw new RuntimeException("Read error: "+path);}
		if(version!=VERSION){fail("Wrong seed-offset version: expected "+VERSION+", found "+version);}
		if(k!=expectedK){fail("Wrong seed K: expected "+expectedK+", found "+k);}
		if(trainingGenes<=0){fail("Missing or invalid training-gene count; require >0, found "+trainingGenes);}
		if(!validSha80(seedSetSha80)){fail("Missing or invalid seed-set sha80: "+seedSetSha80);}
		if(expectedSeedSetSha80!=null && !expectedSeedSetSha80.equals(seedSetSha80)){
			fail("Seed-set sha80 mismatch: expected "+expectedSeedSetSha80+", found "+seedSetSha80);
		}
		if(!sawHeader){fail("Missing exact data header: "+HEADER);}
		for(long seed : expectedSeeds.toArray()){
			if(!map.contains(seed)){fail("Missing expected seed row: "+seed);}
		}
		if(declaredSeeds!=expectedSeeds.size()){
			fail("Declared seed count mismatch: expected "+expectedSeeds.size()+", found "+declaredSeeds);
		}
		if(map.size()!=expectedSeeds.size()){
			fail("Seed-row count mismatch: expected "+expectedSeeds.size()+", found "+map.size());
		}
		return new SeedOffsetTable(k, trainingGenes, seedSetSha80, map);
	}

	private SeedOffsetTable(int k_, long trainingGenes_, String seedSetSha80_, LongObjectMap<KmerInfo> map_){
		k=k_; trainingGenes=trainingGenes_; seedSetSha80=seedSetSha80_; map=map_;
	}

	public KmerInfo get(long seed){return map.get(seed);}
	public int size(){return map.size();}
	public int k(){return k;}
	public long trainingGenes(){return trainingGenes;}
	public String seedSetSha80(){return seedSetSha80;}

	static boolean validSha80(String value){return value!=null && value.matches("[0-9a-f]{20}");}
	private static String ascii(byte[] bytes){return new String(bytes, StandardCharsets.US_ASCII);}
	private static void fail(String message){throw new IllegalArgumentException("Invalid seed-offset table: "+message);}

	/** One immutable resource record; count zero means the seed casts no vote. */
	public static final class KmerInfo {
		KmerInfo(double leftOffset_, double rightOffset_, double leftSD_, double rightSD_, long count_){
			leftOffset=leftOffset_; rightOffset=rightOffset_; leftSD=leftSD_; rightSD=rightSD_; count=count_;
		}
		public boolean trained(){return count>0;}
		public final double leftOffset, rightOffset, leftSD, rightSD;
		public final long count;
	}

	static final int VERSION=1;
	static final String HEADER="#Sequence\tLeftOffset\tRightOffset\tLeftSD\tRightSD\tCount";
	private static final byte[] HEADER_BYTES=HEADER.getBytes(StandardCharsets.US_ASCII);
	private final int k;
	private final long trainingGenes;
	private final String seedSetSha80;
	private final LongObjectMap<KmerInfo> map;
}
