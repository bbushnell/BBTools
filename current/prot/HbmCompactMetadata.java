package prot;

import java.util.HashMap;
import java.util.HashSet;

import fileIO.ByteFile;
import parse.LineParser1;
import structures.ByteBuilder;

/** Optional per-family compact-HBM annotations. Unknown values are never emitted.
 * Identity uses percent, BLOSUM uses the raw alignment score, and HBM uses the
 * path-relative score. Coordinates are zero-based and inclusive on the consensus.
 * @author Collei
 */
public final class HbmCompactMetadata {

	/** Sentinel values are field-specific; negative real score cutoffs remain valid. */
	public double identity=-1, blosum=-999999, hbm=-999999;
	public int start=-1, stop=-1, minlen=-1, maxlen=-1;

	/** Named fields allow a partial cutoff row without writing placeholder values. */
	void append(ByteBuilder out){
		validate();
		if(identity!=-1 || blosum!=-999999 || hbm!=-999999){
			out.append("#cutoff");
			if(identity!=-1){out.append("\tidentity=").append(Double.toString(identity));}
			if(blosum!=-999999){out.append("\tblosum=").append(Double.toString(blosum));}
			if(hbm!=-999999){out.append("\thbm=").append(Double.toString(hbm));}
			out.nl();
		}
		append(out, "#start", start); append(out, "#stop", stop);
		append(out, "#minlen", minlen); append(out, "#maxlen", maxlen);
	}

	private static void append(ByteBuilder out, String name, int value){
		if(value!=-1){out.append(name).tab().append(value).nl();}
	}

	void parse(LineParser1 row){
		final String key=row.parseString(0);
		if(!seen.add(key)){throw bad("duplicate "+key);}
		if(key.equals("#cutoff")){
			if(row.terms()<2 || row.terms()>4){throw bad("empty or oversized cutoff row");}
			final HashSet<String> fields=new HashSet<String>();
			for(int i=1; i<row.terms(); i++){
				final String term=row.parseString(i); final int eq=term.indexOf('=');
				if(eq<1){throw bad("cutoff requires named identity/blosum/hbm fields");}
				final String name=term.substring(0, eq);
				if(!fields.add(name)){throw bad("duplicate cutoff "+name);}
				final double value=Double.parseDouble(term.substring(eq+1));
				if(name.equals("identity") && value!=-1){identity=value;}
				else if(name.equals("blosum") && value!=-999999){blosum=value;}
				else if(name.equals("hbm") && value!=-999999){hbm=value;}
				else{throw bad("unknown cutoff or serialized placeholder: "+term);}
			}
		}else{
			if(row.terms()!=2){throw bad("expected one value for "+key);}
			final int value=row.parseInt(1);
			if(value<0){throw bad("serialized placeholder: "+key);}
			if(key.equals("#start")){start=value;}
			else if(key.equals("#stop")){stop=value;}
			else if(key.equals("#minlen")){minlen=value;}
			else if(key.equals("#maxlen")){maxlen=value;}
			else{throw bad("unknown optional header "+key);}
		}
		validate();
	}

	void validate(){
		if(!Double.isFinite(identity) || !Double.isFinite(blosum) || !Double.isFinite(hbm)
				|| (identity!=-1 && (identity<0 || identity>100))
				|| start< -1 || stop< -1 || minlen< -1 || maxlen< -1
				|| (start>=0 && stop>=0 && start>stop)
				|| minlen==0 || maxlen==0 || (minlen>0 && maxlen>0 && minlen>maxlen)){
			throw bad("invalid cutoff, coordinate or length bounds");
		}
	}

	/** Imports only the established effective gates; endpoint percentiles are absent. */
	static HashMap<String,HbmCompactMetadata> fromProfile(String path){
		final HashMap<String,HbmCompactMetadata> result=new HashMap<String,HbmCompactMetadata>();
		final ByteFile input=ByteFile.makeByteFile(path, false); final LineParser1 row=new LineParser1('\t');
		boolean header=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length>0 && line[0]=='#'){continue;}
				row.set(line);
				if(!header){
					if(row.terms()!=33 || !row.termEquals("rep_id", 2) || !row.termEquals("raw_score_effective", 7)
							|| !row.termEquals("identity_effective", 10) || !row.termEquals("length_lo_effective", 19)
							|| !row.termEquals("length_hi_effective", 22) || !row.termEquals("hbm_path_effective", 31)){
						throw bad("expected schema7 effective threshold columns");
					}
					header=true; continue;
				}
				if(row.terms()!=33 || row.parseInt(0)!=result.size()){throw bad("profile row width/rank mismatch");}
				final HbmCompactMetadata metadata=new HbmCompactMetadata();
				// Schema7FamilyAcceptanceLoader uses canonicalFloat for these two gates.
				// Widen the parsed float exactly; direct decimal-to-double can move a boundary.
				metadata.identity=row.parseFloat(10); metadata.blosum=row.parseInt(7);
				metadata.hbm=row.parseFloat(31); metadata.minlen=row.parseInt(19); metadata.maxlen=row.parseInt(22);
				metadata.validate();
				if(result.put(row.parseString(2), metadata)!=null){throw bad("duplicate profile representative");}
			}
		}finally{if(input.close()){throw bad("profile I/O failure");}}
		if(!header || result.isEmpty()){throw bad("empty profile");}
		return result;
	}

	private static IllegalArgumentException bad(String reason){return new IllegalArgumentException("Compact HBM metadata: "+reason);}
	private final HashSet<String> seen=new HashSet<String>();
}
