package org.panda.resource;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Base class for resource accessor classes that download text files (can be compressed), and loads themselves reading
 * these files.
 *
 * @author Ozgun Babur
 */
public abstract class FileServer
{
	public abstract String[] getLocalFilenames();

	/**
	 * This default method assumes that the distant file is in the GitHub repo and has the same name. If that is
	 * violated, then this method should be overridden.
	 */
	public String[] getDistantURLs()
	{
		return Arrays.stream(getLocalFilenames()).map(name -> GITHUB_REPO_BASE + name).toArray(String[]::new);
	}

	public static final String GITHUB_REPO_BASE =
		"https://raw.githubusercontent.com/PathwayAndDataAnalysis/repo/master/resource-files/";

	public FileServer()
	{
		try
		{
			if (!init()) throw new RuntimeException("Cannot initialize the resource");
		}
		catch (IOException e)
		{
			throw new RuntimeException("Cannot initialize the resource.", e);
		}
	}

	public synchronized boolean init() throws IOException
	{
		if (!localResourceExists())
		{
			if (!downloadResources()) return false;
			if (!processTheDownloadedFiles()) return false;
		}

		return load();
	}

	public boolean downloadResources()
	{
		String loc = ResourceDirectory.get();
		String[] url = getDistantURLs();
		String[] filename = getLocalFilenames();

		for (int i = 0; i < url.length; i++)
		{
			boolean success;

			if (url[i].endsWith(".zip") || url[i].endsWith(".gz"))
			{
				success = Download.downloadAndUncompress(url[i], loc + File.separator + filename[i]);
			}
			else
			{
				success = Download.downloadAsIs(url[i], loc + File.separator + filename[i]);
			}

			if (!success)
			{
				System.err.println("Cannot download " + url[i]);
				return false;
			}
		}
		return true;
	}

	/**
	 * THis method is a hook point for the resource accessor that needs to modify the downloaded file before making it
	 * ready for loading.
	 * @return
	 */
	public boolean processTheDownloadedFiles()
	{
		return true;
	}

	public abstract boolean load() throws IOException;

	public Stream<String> getResourceAsStream(String filename) throws IOException
	{
		return Files.lines(Paths.get(ResourceDirectory.get() + File.separator + filename));
	}

	/**
	 * Override this method if checking id the resource is ready involves more than
	 * @return
	 */
	public boolean localResourceExists()
	{
		for (String filename : getLocalFilenames())
		{
			if (!(new File(ResourceDirectory.get() + File.separator + filename).exists()))
				return false;
		}
		return true;
	}

	public String locateInBase(String file)
	{
		return ResourceDirectory.get() + File.separator + file;
	}
}
