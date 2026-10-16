package prok;

import fileIO.ByteFile;
import map.ObjectMap;
import parse.LineParser1;

/** Versioned loose-resource bindings. Digests describe shipped bytes; callers validate
 * biological bindings, while release fixtures verify the payload digests.
 * @author Brian, Raiden */
final class FamilyResourceManifest {

	FamilyResourceManifest(String path, NcrnaFamilyRegistry.Entry family){
		assert(path!=null && family!=null) : "A manifest must bind one explicit registered family";
		final ByteFile in=RrnaResourceIO.open(path);
		final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "#ncrna_family\t1\t"+family.id, path);
			RrnaResourceIO.header(in, "role\tmodel\tend\tpath\tformat\tbytes\tsha80", path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==7, "Family manifest requires seven columns: "+path);
				final String role=p.parseString(0), model=p.parseString(1), end=p.parseString(2);
				final String relative=p.parseString(3), format=p.parseString(4), digest=p.parseString(6);
				final long bytes=p.parseLong(5);
				require(bytes>0 && digest.matches("[0-9a-f]{20}"), "Resource provenance requires positive bytes and sha80: "+relative);
				validateRelative(relative);
				if(role.equals("table")){
					boolean known=false;for(String name:family.models){known|=name.equals(model);}
					require(known && (end.equals("5prime") || end.equals("3prime")), "Table binding must name a registered model/end: "+model+"/"+end);
				}else{
					require(model.equals(".") && end.equals("."), "Only table roles have model/end bindings");
					require(role.equals("consensus") || role.equals("models") || role.equals("seeds") || role.equals("tables") || role.equals("endpoints"),
						"Unknown resource role for "+family.id+": "+role);
				}
				final String key=key(role, model, end);
				require(entries.get(key)==null, "Duplicate family resource role: "+key);
				require(!paths.contains(relative), "Two roles cannot silently alias one pilot payload: "+relative);
				paths.put(relative, 1);
				validateFormat(role, relative, format);
				entries.put(key, new Entry(relative, format, bytes, digest));
			}
		}finally{require(!in.close(), "Family manifest read failed: "+path);}
	}
	Entry get(String role, String model, String end){
		final Entry e=entries.get(key(role, model, end));
		require(e!=null, "Missing required family manifest role: "+key(role, model, end));return e;
	}
	static String key(String role, String model, String end){return role+":"+model+":"+end;}
	static void validateRelative(String path){
		require(path!=null && !path.isEmpty() && !path.startsWith("/") && path.indexOf('\\')<0 && path.indexOf(':')<0,
			"Family entries must be relative paths: "+path);
		for(String part:path.split("/", -1)){require(!part.isEmpty() && !part.equals(".") && !part.equals(".."), "Family entry escapes its directory: "+path);}
	}
	private static void validateFormat(String role, String path, String format){
		final String expected=role.equals("consensus") || role.equals("seeds") ? "fasta" : role.equals("models") ? "hbm"
			: role.equals("endpoints") ? "bbnets" : "tsv";
		require(format.equals(expected), "Unsupported format for "+role+": "+format);
		if(role.equals("models")){require(path.endsWith(".hbm") || path.endsWith(".hbmc") || path.endsWith(".hbmc.gz"),
			"HBM format dispatch depends on the native .hbm/.hbmc[.gz] suffix: "+path);}
	}
	static final class Entry {
		Entry(String path_, String format_, long bytes_, String digest_){path=path_;format=format_;bytes=bytes_;sha80=digest_;}
		final String path, format, sha80;
		final long bytes;
	}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final ObjectMap<String, Entry> entries=new ObjectMap<String, Entry>(String.class, Entry.class);
	private final map.ObjectIntMap<String> paths=new map.ObjectIntMap<String>(String.class);
}
