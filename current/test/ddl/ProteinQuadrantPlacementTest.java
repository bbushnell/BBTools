package test.ddl;

import java.io.DataOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

import ddl.DDLFamilyProfiles;
import ddl.DDLFamilySearch;
import ddl.DDLKeyTieBreak;
import ddl.DDLQueryOrigins;
import ddl.FamilySketchBuilder;
import ddl.ProteinKmerGeometry;
import ddl.ProteinPositionMode;
import dna.AminoAcid;
import prot.Blosum62;
import prot.cluster.ClusterSketchCodec;
import shared.Tools;

/** Independent word-by-word quadrant oracle and parent-compatible default-byte capture.
 * Run the same fixture against the immutable parent and candidate, then compare default.bin.
 * @author Brian Bushnell, Keqing
 */
public final class ProteinQuadrantPlacementTest{

	public static void main(String[] args)throws Exception{
		check(args.length==2 && (args[0].equals("baseline") || args[0].equals("candidate")),
			"Require baseline|candidate and a fresh fixture output directory");
		final Path out=Paths.get(args[1]);Files.createDirectory(out);
		captureDefaults(out);
		if(args[0].equals("candidate")){
			final ProteinPositionMode quadrant=ProteinPositionMode.parse("fixed-quadrant");
			check(ProteinPositionMode.fromSketch("positional-fixed-quadrant-v1")==quadrant && quadrant!=ProteinPositionMode.FIXED,
				"Quadrant wire identity must not alias the accepted fixed layout");
			metadataIdentity(quadrant);
			for(int width:new int[]{512,1024}){
				for(int overlap:new int[]{0,2,8}){
					for(int length:new int[]{0,4,5,40,49,50,53,54,98,100,151,200,231}){
						compare(sequence(length),width,overlap,quadrant);
					}
					for(int distance:new int[]{0,47,48,49,50,97,98,99,146,147,148,195}){
						final byte[] isolated=new byte[200];Arrays.fill(isolated,(byte)'X');
						System.arraycopy("ACDEF".getBytes(StandardCharsets.US_ASCII),0,isolated,195-distance,5);
						compare(isolated,width,overlap,quadrant);
					}
				}
			}
			for(int partitions:new int[]{2,8,16}){
				boolean rejected=false;
				try{DDLFamilyProfiles.sketchPositional(new byte[0],K5,0,partitions,0,quadrant,new char[1024],null);}
				catch(IllegalArgumentException expected){rejected=true;}
				check(rejected,"Even empty queries must reject unsupported quadrant partition counts");
			}
			boolean rejected=false;
			try{DDLFamilyProfiles.sketchPositional(new byte[0],K5,0,4,0,quadrant,new char[256],null);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected,"Quadrant placement must reject unsupported register widths");
			check(changed>0 && duplicated>0,"Oracle cases must exercise changed buckets and duplicated border words");
		}
		Files.write(out.resolve("PASS"),("PROTEIN_QUADRANT_PLACEMENT_PASS phase="+args[0]+" checks="+checks+
			" changed_bins="+changed+" border_copies="+duplicated+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Distinct wire identity must survive headers and reject missing or cross-mode bindings. */
	private static void metadataIdentity(ProteinPositionMode mode){
		final HashMap<String,String> meta=new HashMap<String,String>();
		DDLFamilySearch.validatePositionalMetadata(meta,ProteinPositionMode.FIXED.sketch,4,8,K5);
		reject(()->DDLFamilySearch.validatePositionalMetadata(meta,mode.sketch,4,8,K5),
			"Only legacy FIXED may infer missing four/eight metadata");
		meta.put("partitions","4");meta.put("overlap","2");meta.put("tag_bits","2");meta.put("segment_words","49");
		DDLFamilySearch.validatePositionalMetadata(meta,mode.sketch,4,2,K5);
		meta.put("position_basis","residue-length");
		reject(()->DDLFamilySearch.validatePositionalMetadata(meta,mode.sketch,4,2,K5),
			"Quadrant placement retains fixed coordinates, never relative metadata");
		meta.remove("position_basis");meta.put("partitions","8");
		reject(()->DDLFamilySearch.validatePositionalMetadata(meta,mode.sketch,8,2,K5),
			"Stored quadrant identity cannot accept an eight-partition layout");
		final ClusterSketchCodec.Header sliced=new ClusterSketchCodec.Header(K5,512,73,mode,4,2,1,null);
		final ClusterSketchCodec.Header fixed=new ClusterSketchCodec.Header(K5,512,73,ProteinPositionMode.FIXED,4,2,1,null);
		check(!sliced.sameGeometry(fixed) && !fixed.sameGeometry(sliced) && sliced.withFamilies(2).positionMode==mode,
			"Codec geometry must distinguish placement and preserve its identity when roster size changes");
		final DDLKeyTieBreak.PolicyContext policy=new DDLKeyTieBreak.PolicyContext(DDLKeyTieBreak.Policy.SMALLEST,null,
			512,K5,73,mode,4,2,null,null);
		reject(()->policy.requireGeometry(512,K5,73,ProteinPositionMode.FIXED,4,2,null),
			"A selection context cannot be reused across incompatible bucket placement");
		reject(()->ProteinPositionMode.resolve("fixed-quadrant","positional-v1",null,ProteinPositionMode.PLAIN),
			"Contradictory position flags cannot erase quadrant identity");
	}

	/** This method calls only APIs present in the accepted parent and names only old modes. */
	private static void captureDefaults(Path out)throws Exception{
		try(DataOutputStream stream=new DataOutputStream(Files.newOutputStream(out.resolve("default.bin")))){
			for(ProteinPositionMode mode:new ProteinPositionMode[]{ProteinPositionMode.PLAIN,ProteinPositionMode.FIXED,ProteinPositionMode.RELATIVE}){
				for(int width:new int[]{512,1024}){
					for(int overlap:new int[]{0,2,8}){
						final byte[] raw=sequence(200),encoded=Blosum62.encode(raw,"quadrant fixture");
						final FamilySketchBuilder builder=new FamilySketchBuilder(width,K5,73,mode,4,overlap,null);
						final char[] a=new char[width],b=new char[width];
						final long words=builder.sketchRaw(raw,a);
						check(builder.sketchEncoded(encoded,b)==words && Arrays.equals(a,b),"Legacy raw/encoded default bytes differ");
						builder.addEncoded(encoded);builder.writeProfile(b);
						check(Arrays.equals(a,b),"Legacy one-member profile must preserve query bytes");
						stream.writeUTF(mode.sketch);stream.writeInt(width);stream.writeInt(overlap);stream.writeLong(words);
						for(char value:a){stream.writeChar(value);}
					}
				}
			}
		}
	}

	/** Reconstruct each word independently, with literal K5/four-band boundaries; no production placement helper. */
	private static void compare(byte[] raw,int width,int overlap,ProteinPositionMode mode){
		final char[] expected=new char[width],actual=new char[width],encoded=new char[width],traced=new char[width],legacy=new char[width];
		final int[] first=new int[width],last=new int[width],tags=new int[width];Arrays.fill(first,-1);Arrays.fill(last,-1);
		final DDLKeyTieBreak.KeySet keys=new DDLKeyTieBreak.KeySet(width),expectedKeys=new DDLKeyTieBreak.KeySet(width),encodedKeys=new DDLKeyTieBreak.KeySet(width);
		final byte[] codes=Blosum62.encode(raw,"quadrant fixture");
		final long enumerated=DDLFamilyProfiles.enumerateKeys(raw,K5,73,mode,4,overlap,null,keys);
		final long encodedWords=DDLFamilyProfiles.enumerateKeysEncoded(codes,K5,73,mode,4,overlap,null,encodedKeys);
		long words=0;
		for(int start=0;start+5<=raw.length;start++){
			long word=0;boolean valid=true;
			for(int j=0;j<5;j++){
				final int aa=AminoAcid.acidToNumberNoStops[raw[start+j]];
				if(aa<0){valid=false;break;}
				word=(word<<5)|aa;
			}
			if(!valid){continue;}words++;
			final int distance=raw.length-5-start,primary=Math.min(3,distance/49),maxTag=Math.min(3,(raw.length-5)/49);
			for(int tag=0;tag<=maxTag;tag++){
				if(tag!=primary && !(tag==primary-1 && distance<primary*49+overlap) &&
					!(tag==primary+1 && distance>=(primary+1)*49-overlap)){continue;}
				if(tag!=primary){duplicated++;}
				final long hash=Tools.hash64shift(((word<<2)|tag)^73L);
				final int zeros=Long.numberOfLeadingZeros(hash)&63;
				final int key=Math.max(1,(zeros<<10)+(int)((~(hash>>>(53-zeros)))&1023));
				final int bucket=tag*(width/4)+(int)(hash%(width/4)+(width/4))%(width/4);
				check(keys.contains(bucket,key),"All-key enumeration lost the original tagged key in its required quadrant");
				check(encodedKeys.contains(bucket,key),"Encoded all-key enumeration lost an oracle key");
				expectedKeys.add(bucket,key);
				if(key>expected[bucket]){expected[bucket]=(char)key;first[bucket]=last[bucket]=start;tags[bucket]=1<<tag;}
				else if(key==expected[bucket]){last[bucket]=start;tags[bucket]|=1<<tag;}
			}
		}
		check(enumerated==words && keys.size()==expectedKeys.size(),"All-key enumeration has missing or extra keys/word counts");
		final FamilySketchBuilder builder=new FamilySketchBuilder(width,K5,73,mode,4,overlap,null);
		check(builder.sketchRaw(raw,actual)==words && Arrays.equals(actual,expected),"Quadrant sketch differs from independent placement oracle");
		check(builder.sketchEncoded(codes,encoded)==words && Arrays.equals(actual,encoded),"Encoded member changed quadrant placement");
		check(encodedWords==words && encodedKeys.size()==expectedKeys.size(),
			"Encoded all-key enumeration changed counts");
		encodedKeys.maxima(encoded);check(Arrays.equals(encoded,expected),"Encoded enumeration maxima differ from oracle");
		final DDLQueryOrigins origins=new DDLQueryOrigins(width);
		check(DDLFamilyProfiles.sketchPositional(raw,K5,73,4,overlap,traced,null,origins,mode)==words && Arrays.equals(traced,expected),
			"Origin tracing changed keys or word counts");
		check(Arrays.equals(first,origins.first) && Arrays.equals(last,origins.last) && Arrays.equals(tags,origins.tags),
			"Origin coordinates or quadrant tags differ from independent winners");
		builder.addEncoded(codes);builder.writeProfile(encoded);check(Arrays.equals(encoded,expected),"Family member lost quadrant placement");
		check(builder.fork().positionMode()==mode && !builder.matchesGeometry(width,K5,73,ProteinPositionMode.FIXED,4,overlap,null),
			"Worker fork or geometry identity can mix legacy and quadrant layouts");
		DDLFamilyProfiles.sketchPositional(raw,K5,73,4,overlap,ProteinPositionMode.FIXED,legacy,null);
		for(int i=0;i<width;i++){if(legacy[i]!=actual[i]){changed++;}}
		check(DDLFamilyProfiles.positionMetadataColumns(mode).equals(DDLFamilyProfiles.positionMetadataColumns(ProteinPositionMode.FIXED)) &&
			DDLFamilyProfiles.positionMetadataValues(mode,4,overlap,K5).equals(DDLFamilyProfiles.positionMetadataValues(ProteinPositionMode.FIXED,4,overlap,K5)),
			"Quadrant placement must retain fixed coordinate metadata");
	}

	private static byte[] sequence(int length){
		final byte[] out=new byte[length],letters="ACDEFGHIKLMNPQRSTVWY".getBytes(StandardCharsets.US_ASCII);
		for(int i=0;i<length;i++){out[i]=i%31==17 ? (byte)'X' : letters[(i*i+7*i+3)%letters.length];}
		return out;
	}
	private static void check(boolean value,String message){checks++;if(!value){throw new AssertionError(message);}}
	private static void reject(Runnable action,String message){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed,message);
	}
	private static final ProteinKmerGeometry K5=ProteinKmerGeometry.exact(5);
	private static long checks,changed,duplicated;
}
