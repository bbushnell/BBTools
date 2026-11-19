package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.nio.charset.StandardCharsets;

import parse.LineParser1;

/** Reuses an explicitly supplied prior digest only for a sealed, unchanged composite.
 * The producer must derive this small pinned record from accepted verification evidence;
 * recording current metadata beside an unverified digest does not establish that evidence.
 * @author Brian Bushnell, Yoimiya */
public final class MagQCVerifiedNetMetadata {

	private MagQCVerifiedNetMetadata(){}

	/** No payload read; fail closed on any identity, timestamp, permission or binding change. */
	public static void verify(final String manifest, final String manifestPin, final String network, final String networkPin){
		try{
			final Path metadata=Paths.get(manifest), path=Paths.get(network).toAbsolutePath().normalize();
			require(Files.isRegularFile(metadata) && Files.size(metadata)<=65536,
				"Verified-net metadata must be a small regular record, never a model payload");
			final byte[] content;
			try(java.io.InputStream in=Files.newInputStream(metadata)){
				final byte[] buffer=new byte[65537];int used=0,n;
				while(used<buffer.length && (n=in.read(buffer,used,buffer.length-used))>=0){used+=n;}
				require(used<=65536,"Verified-net metadata grew beyond its bounded record size");content=Arrays.copyOf(buffer,used);
			}
			require(DigestSuffix.bytes(content).equals(DigestSuffix.requireSuffix(manifestPin,"net metadata pin")),
				"Verified-net metadata differs from its independently supplied pin");
			final LineParser1 p=new LineParser1('\t');
			int first=-1,last=-1;for(int i=0;i<content.length;i++){if(content[i]=='\n'){if(first<0){first=i;}else{require(last<0,"Verified-net metadata has extra rows");last=i;}}}
			require(first>=0 && last==content.length-1 && last>first+1,"Verified-net metadata needs a header and one complete row");
			{
				require(Arrays.equals(Arrays.copyOfRange(content,0,first),HEADER.getBytes(StandardCharsets.US_ASCII)),"Unsupported verified-net metadata header");
				p.set(Arrays.copyOfRange(content,first+1,last));
				require(p.terms()==7,"Verified-net metadata must contain exactly one seven-field row");
				final Path recorded=Paths.get(p.parseString(0));
				require(recorded.isAbsolute() && recorded.normalize().equals(path),"Verified-net path differs from the configured network");
				require(p.termEquals(DigestSuffix.requireSuffix(networkPin,"configured net pin"),1),"Verified-net digest differs from the configured network pin");
				require(!Files.isSymbolicLink(path),"Verified-net metadata does not authorize a replacement file symlink");
				final BasicFileAttributes a=Files.readAttributes(path,BasicFileAttributes.class);
				final FileTime ctime=(FileTime)Files.getAttribute(path,"unix:ctime");
				final long inode=((Number)Files.getAttribute(path,"unix:ino")).longValue();
				final int mode=((Number)Files.getAttribute(path,"unix:mode")).intValue()&0777;
				require(a.isRegularFile() && a.size()==p.parseLong(2),"Verified-net size or regular-file identity changed");
				require(a.lastModifiedTime().to(TimeUnit.NANOSECONDS)==p.parseLong(3)
					&& ctime.to(TimeUnit.NANOSECONDS)==p.parseLong(4),"Verified-net modification or change time changed");
				// Linux inode values are unsigned; device numbers vary across Dori mounts.
				require(Long.toUnsignedString(inode).equals(p.parseString(5)),"Verified-net inode changed");
				require(mode==Integer.parseInt(p.parseString(6),8) && (mode&0222)==0,"Verified-net permissions changed or permit writing");
			}
		}catch(java.io.IOException e){throw new IllegalArgumentException("Cannot verify accepted network metadata",e);}
	}

	private static void require(final boolean ok,final String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	public static final String HEADER="path\tsha80\tbytes\tmtime_ns\tctime_ns\tinode\tmode";
}
