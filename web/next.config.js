/** @type {import('next').NextConfig} */
const nextConfig = {
  images: { unoptimized: true },
  trailingSlash: false,
  // Vercel: remove output: 'export' — Vercel handles server-side rendering
  // Static export: add `output: 'export'` and `trailingSlash: true` for `npx next build && npx serve out`
}

module.exports = nextConfig