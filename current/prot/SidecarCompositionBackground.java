package prot;

import java.util.HashMap;

import fileIO.ByteFile;
import parse.LineParser1;

/** Reads a pinned sidecar's composition background without loading its family index.
 * Keeps the serialized vector intact so independent packs use exactly the same
 * double values after loading. This does not validate the donor's family rows.
 * @author Keqing
 */
final class SidecarCompositionBackground {

	static SidecarCompositionBackground load(final String path, final String pin,
			final String variant, final double pseudocount, final String mode){
		DigestSuffix.requireSuffix(pin, "backgroundsha80");
		require(Double.isFinite(pseudocount) && pseudocount>=0, "Inherited composition pseudocount must be finite and nonnegative");
		require(pin.equals(DigestSuffix.file(path)), "Source background file does not match backgroundsha80");
		final HashMap<String,String> header=new HashMap<String,String>();
		final ByteFile bf=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		byte[] vectorLine=null, pseudocountLine=null;
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0){continue;}
				if(line[0]!='#'){break;}
				lp.set(line);
				final String key=lp.parseString(0);
				if(key.startsWith("#background_") || key.startsWith("#composition_")
						|| key.equals("#kind") || key.equals("#schema_version")){
					require(lp.terms()==2, "Background metadata must have exactly two columns: "+key);
					require(!key.startsWith("#background_") || key.equals("#background_400f"), "Unsupported source background vector: "+key);
					require(header.put(key, lp.parseString(1))==null, "Duplicate background metadata: "+key);
				}
				if(key.equals("#background_400f")){vectorLine=line.clone();}
				if(key.equals("#composition_rawpseudo")){pseudocountLine=line.clone();}
			}
		}finally{
			require(!bf.close(), "I/O error reading source background: "+path);
		}
		require("family_shortlist_sidecar_v1".equals(header.get("#kind")), "Source must be a family shortlist sidecar");
		final String schema=header.get("#schema_version");
		require("1".equals(schema) || "2".equals(schema), "Unsupported source sidecar schema");
		require("aa20".equals(header.get("#composition_alphabet"))
			&& "2".equals(header.get("#composition_k"))
			&& "400".equals(header.get("#composition_dims")), "Source must use aa20 dimers with 400 dimensions");
		require(variant.equals(header.get("#composition_variant"))
			&& mode.equals(header.get("#composition_background")), "Source composition transform/background mode differs from requested recipe");
		final String raw=header.get("#composition_rawpseudo");
		require(raw!=null && Double.parseDouble(raw)==pseudocount, "Source composition pseudocount differs from requested recipe");
		require(vectorLine!=null, "Source sidecar lacks its 400-dimensional background vector");
		lp.set(vectorLine);
		final LineParser1 values=new LineParser1(',');
		values.set(lp.parseByteArray(1));
		require(values.terms()==400, "Source background must have exactly 400 probabilities");
		final double[] frequencies=new double[400];
		double sum=0;
		for(int i=0; i<frequencies.length; i++){
			// Match FamilyShortlistSidecar.load's decoder exactly. This is a
			// one-time 400-value header, not a per-family or per-query hot path.
			final double value=Double.parseDouble(values.parseString(i));
			require(Double.isFinite(value) && value>=0 && value<=1, "Invalid source background probability at dimension "+i);
			frequencies[i]=value; sum+=value;
		}
		// Existing sidecars round each component to eight decimal places; do not
		// renormalize them because the consumer compares their loaded values exactly.
		require(Math.abs(sum-1)<=400*5.1e-9, "Source background probabilities do not sum to one within serialization precision");
		return new SidecarCompositionBackground(frequencies, vectorLine, pseudocountLine, pin);
	}

	private SidecarCompositionBackground(final double[] frequencies_, final byte[] vectorLine_,
			final byte[] pseudocountLine_, final String pin_){
		assert(frequencies_.length==400 && vectorLine_!=null && pseudocountLine_!=null) : "Builder and consumer require the same serialized background and pseudocount";
		frequencies=frequencies_; vectorLine=vectorLine_; pseudocountLine=pseudocountLine_; pin=pin_;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	final double[] frequencies;
	final byte[] vectorLine, pseudocountLine;
	final String pin;
}
