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
		final NcrnaFamilyRegistry.Group group=family.group;
		final int shift=group==null ? 0 : 1;
		try{
			RrnaResourceIO.header(in, group==null ? "#ncrna_family\t1\t"+family.id : "#ncrna_family\t2\t"+group.id, path);
			RrnaResourceIO.header(in, (group==null ? "" : "member\t")+"role\tmodel\tend\tpath\tformat\tbytes\tsha80", path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==7+shift, (group==null ? "Family manifest requires seven columns: " : "Group manifest requires eight columns: ")+path);
				final String member=group==null ? "" : p.parseString(0);
				final String role=p.parseString(shift), model=p.parseString(shift+1), end=p.parseString(shift+2);
				final String relative=p.parseString(shift+3), format=p.parseString(shift+4), digest=p.parseString(shift+6);
				final long bytes=p.parseLong(shift+5);
				if(group!=null){
					require(group.hasMember(member) || member.equals(".") && group.sharedRoles.length>0, "Unknown group member "+member+" in "+path);
					require(group.permits(member, role), "Unsupported resource role "+role+" for group "+group.id+" member "+member);
				}
				require(bytes>0 && digest.matches("[0-9a-f]{20}"), "Resource provenance requires positive bytes and sha80: "+relative);
				validateRelative(relative);
				if(role.equals("table")){
					boolean known=false;for(String name:family.models){known|=name.equals(model);}
					require(known && (end.equals("5prime") || end.equals("3prime")), "Table binding must name a registered model/end: "+model+"/"+end);
				}else{
					require(model.equals(".") && end.equals("."), "Only table roles have model/end bindings");
					require(role.equals("consensus") || role.equals("models") || role.equals("seeds") || role.equals("tables") || role.equals("endpoints") || role.equals("cutoffs")
						|| role.equals("boundary") || role.equals("boundary_start") || role.equals("boundary_stop") || role.equals("legacy_consensus")
						|| role.equals("vote_offsets") || role.equals("joined_model_positions") || role.equals("length_bounds"),
						"Unknown resource role for "+family.id+": "+role);
				}
				final String key=key(member, role, model, end);
				require(entries.get(key)==null, "Duplicate family resource role: "+key);
				final String memberPath=member.isEmpty() || member.equals(".") ? relative : member+'/'+relative;
				require(!paths.contains(memberPath), "Two roles cannot silently alias one payload: "+memberPath);
				paths.put(memberPath, 1);
				validateFormat(role, relative, format);
				entries.put(key, new Entry(member, relative, format, bytes, digest));
			}
		}finally{require(!in.close(), "Family manifest read failed: "+path);}
		if(group!=null){
			for(String role:group.sharedRoles){
				require(entries.get(key(".", role, ".", "."))!=null,
					"Missing required shared group resource: "+group.id+"/"+role+" in "+path);
			}
			for(String member:group.members){for(String role:group.requiredRoles){
				require(entries.get(key(member, role, ".", "."))!=null,
					"Missing required group member resource: "+group.id+"/"+member+"/"+role+" in "+path);
			}}
		}
	}
	Entry get(String role, String model, String end){
		return get("", role, model, end);
	}
	Entry get(String member, String role, String model, String end){
		final Entry e=entries.get(key(member, role, model, end));
		require(e!=null, "Missing required family manifest role: "+key(member, role, model, end));return e;
	}
	static String key(String member, String role, String model, String end){return (member.isEmpty() ? "" : member+'/')+key(role, model, end);}
	static String key(String role, String model, String end){return role+":"+model+":"+end;}
	static void validateRelative(String path){
		require(path!=null && !path.isEmpty() && !path.startsWith("/") && path.indexOf('\\')<0 && path.indexOf(':')<0,
			"Family entries must be relative paths: "+path);
		for(String part:path.split("/", -1)){require(!part.isEmpty() && !part.equals(".") && !part.equals(".."), "Family entry escapes its directory: "+path);}
	}
	private static void validateFormat(String role, String path, String format){
		final String expected=role.equals("consensus") || role.equals("legacy_consensus") || role.equals("seeds") ? "fasta" : role.equals("models") ? "hbm"
			: role.equals("endpoints") ? "bbnets" : role.equals("boundary") ? "bbnet" : "tsv";
		require(format.equals(expected), "Unsupported format for "+role+": "+format);
		if(role.equals("models")){require(path.endsWith(".hbm") || path.endsWith(".hbmc") || path.endsWith(".hbmc.gz"),
			"HBM format dispatch depends on the native .hbm/.hbmc[.gz] suffix: "+path);}
	}
	static final class Entry {
		Entry(String member_, String path_, String format_, long bytes_, String digest_){member=member_;path=path_;format=format_;bytes=bytes_;sha80=digest_;}
		final String member, path, format, sha80;
		final long bytes;
	}
	static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final ObjectMap<String, Entry> entries=new ObjectMap<String, Entry>(String.class, Entry.class);
	private final map.ObjectIntMap<String> paths=new map.ObjectIntMap<String>(String.class);
}
