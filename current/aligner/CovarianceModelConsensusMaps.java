package aligner;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.FileFormat;
import map.ObjectIntMap;
import parse.LineParser1;
import structures.ByteBuilder;

/** Loads once-per-consensus placement maps from a saved global CM alignment.
 * Resource/run manifests bind exact CM and Stockholm sha80s. This reader checks
 * names, actual consensus residues, complete RF columns and alignment geometry.
 * Missing rows remain explicit null maps for the caller's exact fallback.
 * @author Brian Bushnell, Raiden
 */
public final class CovarianceModelConsensusMaps {
	public static CovarianceModelConsensusMap[] read(CovarianceModel model, String path,
			String[] names, byte[][] consensus){
		require(model!=null && path!=null && names!=null && consensus!=null && names.length==consensus.length && names.length>0,
			"Placement resources must bind a loaded CM and a nonempty caller consensus library");
		final ObjectIntMap<String> indices=new ObjectIntMap<String>(String.class);
		final ByteBuilder[] rows=new ByteBuilder[names.length];
		for(int i=0; i<names.length; i++){
			final String id=first(names[i]);
			require(consensus[i]!=null && consensus[i].length>0 && !indices.contains(id), "Duplicate name or missing caller consensus: "+id);
			indices.put(id, i);
		}
		final ByteFile in=ByteFile.makeByteFile(FileFormat.testInput(path, FileFormat.TEXT, null, true, true));
		final LineParser1 p=new LineParser1(' ');final ByteBuilder rf=new ByteBuilder();boolean ended=false;
		try{
			require(Arrays.equals(in.nextLine(), "# STOCKHOLM 1.0".getBytes(StandardCharsets.US_ASCII)), "Expected a single global Stockholm alignment");
			for(byte[] line; (line=in.nextLine())!=null;){
				line=compact(line);if(line.length==0){continue;}
				require(!ended, "Unexpected content after Stockholm terminator");
				p.set(line);
				if(p.termEquals("//", 0)){require(p.terms()==1, "Malformed Stockholm terminator");ended=true;continue;}
				if(line[0]=='#'){
					if(p.termEquals("#=GC", 0) && p.terms()>=2 && p.termEquals("RF", 1)){
						require(p.terms()==3, "RF annotation requires one alignment fragment");p.appendTerm(rf, 2);
					}
					continue;
				}
				require(p.terms()==2, "Consensus alignment row requires name and sequence");
				final String id=p.parseString(0);require(indices.contains(id), "Alignment contains a consensus absent from the caller library: "+id);
				final int index=indices.get(id);if(rows[index]==null){rows[index]=new ByteBuilder();}p.appendTerm(rows[index], 1);
			}
		}finally{require(!in.close(), "Consensus alignment read failed: "+path);}
		require(ended && rf.length()>0, "Placement resource requires complete RF annotation and a terminator");
		final byte[] reference=rf.toBytes();final CovarianceModelConsensusMap[] result=new CovarianceModelConsensusMap[names.length];int mapped=0;
		for(int i=0; i<result.length; i++){
			if(rows[i]!=null){result[i]=new CovarianceModelConsensusMap(model, consensus[i], rows[i].toBytes(), reference);mapped++;}
		}
		require(mapped>0, "Empty placement resources cannot silently disable all anchors");return result;
	}
	private static String first(String name){
		require(name!=null && !name.isEmpty(), "Caller consensus needs a stable identifier");
		int end=0;while(end<name.length() && !Character.isWhitespace(name.charAt(end))){end++;}
		require(end>0, "Consensus identifier cannot start with whitespace");return name.substring(0, end);
	}
	private static byte[] compact(byte[] line){
		assert(line!=null):"Whitespace normalization requires a real input row";
		int size=0;boolean space=false;
		for(byte b:line){if(b==' ' || b=='\t' || b=='\r'){space=size>0;}else{if(space){line[size++]=' ';space=false;}line[size++]=b;}}
		return Arrays.copyOf(line, size);
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private CovarianceModelConsensusMaps(){}
}
