/** @type {import('next').NextConfig} */
const isGithubPages = process.env.GITHUB_PAGES === "true"
const repoName = "linear-scaling-pangenome-sv"
const basePath = isGithubPages ? `/${repoName}` : ""

const nextConfig = {
  images: { unoptimized: true },
  trailingSlash: isGithubPages,
  env: {
    NEXT_PUBLIC_BASE_PATH: basePath,
  },
  // Vercel: GITHUB_PAGES is unset, so this block is skipped and Vercel
  // handles server-side rendering normally.
  ...(isGithubPages && {
    output: "export",
    basePath,
    assetPrefix: `${basePath}/`,
  }),
}

module.exports = nextConfig