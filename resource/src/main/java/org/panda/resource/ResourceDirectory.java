package org.panda.resource;

import java.io.File;

/**
 * Many FileServer classes use this class to know which resource directory to use.
 *
 * @author Ozgun Babur
 */
public class ResourceDirectory
{
	public static final String DEFAULT_RESOURCE_DIR_NAME = ".panda";

	private static String directory;

	public static String get()
	{
		if (directory == null)
		{
			String userDir = System.getProperty("user.dir");
			File dir = new File(userDir + File.separator + DEFAULT_RESOURCE_DIR_NAME);
			if (dir.exists() || dir.mkdirs())
			{
				directory = dir.getPath();
			}
			else
			{
				String userHome = System.getProperty("user.home");
				dir = new File(userHome + File.separator + DEFAULT_RESOURCE_DIR_NAME);
				if (dir.exists() || dir.mkdirs()) directory = dir.getPath();
			}
		}
		return directory;
	}

	public static void set(String directory)
	{
		ResourceDirectory.directory = directory;
	}
}
