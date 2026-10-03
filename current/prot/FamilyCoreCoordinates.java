package prot;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;

/** Hash-bound per-family 50%-member-coverage core coordinates for construction gates. */
public final class FamilyCoreCoordinates {

	private static final String DATA_HEADER=
		"rank\trep_id\tconsensus_length\tusable_members\tcoverage_threshold\tcore_start\tcore_end\tcore_length";

	public final int nFamilies;
	public final String[] repIds;
	public final int[] consensusLength,usableMembers,coverageThreshold,coreStart,coreEnd,coreLength;
	public final String familyListSha80,consensusSha80,taskManifestSha80,artifactSha80;

	private FamilyCoreCoordinates(final String[] repIds_, final int[] consensusLength_,
			final int[] usableMembers_, final int[] coverageThreshold_, final int[] coreStart_,
			final int[] coreEnd_, final int[] coreLength_, final String familyListSha80_,
			final String consensusSha80_, final String taskManifestSha80_, final String artifactSha80_){
		repIds=repIds_; nFamilies=repIds.length; consensusLength=consensusLength_;
		usableMembers=usableMembers_; coverageThreshold=coverageThreshold_; coreStart=coreStart_;
		coreEnd=coreEnd_; coreLength=coreLength_; familyListSha80=familyListSha80_;
		consensusSha80=consensusSha80_; taskManifestSha80=taskManifestSha80_; artifactSha80=artifactSha80_;
	}

	public static FamilyCoreCoordinates load(final String path, final String familyListPath,
			final String consensusPath, final String[] expectedRepIds, final int[] expectedConsensusLengths){
		return load(path,familyListPath,consensusPath,expectedRepIds,expectedConsensusLengths,null);
	}

	public static FamilyCoreCoordinates load(final String path, final String familyListPath,
			final String consensusPath, final String[] expectedRepIds, final int[] expectedConsensusLengths,
			final String expectedSemanticFamilySha80){
		if(expectedRepIds==null || expectedConsensusLengths==null || expectedRepIds.length!=expectedConsensusLengths.length){
			throw new IllegalArgumentException("Expected core-coordinate roster/length arrays must be non-null and equal-sized");
		}
		final String before=MagQCTextResource.sha80(path),familyFileSha=MagQCTextResource.sha80(familyListPath),consensusSha=MagQCTextResource.sha80(consensusPath);
		final String expectedFamilySha=expectedSemanticFamilySha80==null?familyFileSha:normalizeSha80(expectedSemanticFamilySha80,"expected semantic family sha80");
		final HashMap<String,String> header=new HashMap<String,String>();
		final ArrayList<String> reps=new ArrayList<String>();
		final ArrayList<Integer> cLen=new ArrayList<Integer>(),usable=new ArrayList<Integer>(),threshold=new ArrayList<Integer>(),
			start=new ArrayList<Integer>(),end=new ArrayList<Integer>(),coreLen=new ArrayList<Integer>();
		final ByteFile bf=ByteFile.makeByteFile(path,false); final LineParser1 lp=new LineParser1('\t');
		boolean sawDataHeader=false; int next=0;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){throw new RuntimeException("Blank core-coordinate line: "+path);}
				if(line[0]=='#'){
					if(sawDataHeader){throw new RuntimeException("Core-coordinate metadata appears after data header");}
					final int tab=indexOf(line,(byte)'\t',1);
					if(tab<2 || tab==line.length-1){throw new RuntimeException("Malformed core-coordinate metadata line");}
					final String key=new String(line,1,tab-1,StandardCharsets.US_ASCII),
						value=new String(line,tab+1,line.length-tab-1,StandardCharsets.US_ASCII);
					if(header.put(key,value)!=null){throw new RuntimeException("Duplicate core-coordinate metadata key: "+key);}
					continue;
				}
				if(!sawDataHeader){
					final String observed=new String(line,StandardCharsets.US_ASCII);
					if(!DATA_HEADER.equals(observed)){throw new RuntimeException("Unexpected core-coordinate data header: "+observed);}
					sawDataHeader=true; continue;
				}
				lp.set(line); if(lp.terms()!=8){throw new RuntimeException("Core-coordinate row "+next+" has "+lp.terms()+" fields, expected 8");}
				final int rank=lp.parseInt(0),cl=lp.parseInt(2),u=lp.parseInt(3),cut=lp.parseInt(4),s=lp.parseInt(5),e=lp.parseInt(6),l=lp.parseInt(7);
				if(rank!=next){throw new RuntimeException("Core-coordinate rank "+rank+" != "+next);}
				if(cl<1 || u<1 || cut!=(u+1)/2 || s<0 || e<s || e>=cl || l!=e-s+1){
					throw new RuntimeException("Invalid core coordinates at rank "+rank+": consensus="+cl+" usable="+u+
						" threshold="+cut+" core=["+s+","+e+"] length="+l);
				}
				reps.add(lp.parseString(1)); cLen.add(Integer.valueOf(cl)); usable.add(Integer.valueOf(u));
				threshold.add(Integer.valueOf(cut)); start.add(Integer.valueOf(s)); end.add(Integer.valueOf(e)); coreLen.add(Integer.valueOf(l)); next++;
			}
		}finally{if(bf.close()){throw new RuntimeException("I/O error reading core coordinates: "+path);}}
		if(!sawDataHeader || next!=expectedRepIds.length){throw new RuntimeException("Core-coordinate family count "+next+" != expected "+expectedRepIds.length);}
		requireHeader(header,"schema_version","1",path); requireHeader(header,"artifact_type","family_core_coordinates",path);
		final String recordedFamily=normalizeSha80(required(header,"family_list_sha80",path),"family_list_sha80");
		final String recordedConsensus=normalizeSha80(required(header,"consensus_sha80",path),"consensus_sha80");
		final String recordedTask=normalizeSha80(required(header,"task_manifest_sha80",path),"task_manifest_sha80");
		if(!recordedFamily.equals(expectedFamilySha)){throw new RuntimeException("Core-coordinate family-list hash mismatch");}
		if(!recordedConsensus.equals(consensusSha)){throw new RuntimeException("Core-coordinate consensus hash mismatch");}
		final String[] repArray=reps.toArray(new String[reps.size()]);
		for(int i=0; i<repArray.length; i++){
			if(!repArray[i].equals(expectedRepIds[i])){throw new RuntimeException("Core-coordinate rep_id mismatch at rank "+i);}
			if(cLen.get(i).intValue()!=expectedConsensusLengths[i]){throw new RuntimeException("Core-coordinate consensus length mismatch at rank "+i);}
		}
		if(!before.equals(MagQCTextResource.sha80(path)) || !familyFileSha.equals(MagQCTextResource.sha80(familyListPath)) ||
				!consensusSha.equals(MagQCTextResource.sha80(consensusPath))){throw new RuntimeException("Core-coordinate resource changed during load");}
		return new FamilyCoreCoordinates(repArray,toInts(cLen),toInts(usable),toInts(threshold),toInts(start),toInts(end),toInts(coreLen),recordedFamily,recordedConsensus,recordedTask,before);
	}

	private static int[] toInts(final ArrayList<Integer> list){final int[] out=new int[list.size()];for(int i=0;i<out.length;i++){out[i]=list.get(i).intValue();}return out;}
	private static int indexOf(final byte[] a,final byte b,final int from){for(int i=from;i<a.length;i++){if(a[i]==b){return i;}}return -1;}
	private static String required(final HashMap<String,String> map,final String key,final String path){final String value=map.get(key);if(value==null||value.isEmpty()){throw new RuntimeException("Missing #"+key+" in "+path);}return value;}
	private static void requireHeader(final HashMap<String,String> map,final String key,final String expected,final String path){final String value=required(map,key,path);if(!value.equals(expected)){throw new RuntimeException("#"+key+"="+value+" != "+expected+" in "+path);}}
	private static String normalizeSha80(final String value,final String label){if(value.length()!=DigestSuffix.HEX_LENGTH){throw new RuntimeException(label+" must be sha80");}for(int i=0;i<value.length();i++){final char c=value.charAt(i);if(!((c>='0'&&c<='9')||(c>='a'&&c<='f'))){throw new RuntimeException(label+" is not lowercase sha80");}}return value;}
}
